#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include<chrono>
#include<string>
#include"mesh.hpp"
#include"state.hpp"
#include"roe.hpp"
#include"roe_I.hpp"
#include"roe_II.hpp"
#include"utilities.hpp"
#include"time_solver.hpp"

using namespace std;

int main()
{
	/********* Domain definition ***********/
	double xa(0), xb(1), dx(1e-2), T(100), t(0), cfl(0.5);
	mesh M (xa, xb, dx);
	string path = "results/"; //"../../results/Euler_3x3_q1d/err_extrapol/isentropic/diff_ord1/dx5e-3/big_da/da005/";
	ofstream fout_param(path+"param_istr.dat");//, ios::out | ios::app);
	ofstream fout_J(path+"J_istr.dat");//, ios::out | ios::app);
	ofstream fout_gradJ(path+"gradJ_istr.dat");
	int N = M.get_N();
	/***************************************/
	
	/************* Shock tube **************/
/*	double uL(0.75), uR(0), rhoL(1), rhoR(0.125), pL(1.0), pR(0.1), x_c(0.5), gamma(1.4);
	double s_uL(0), s_uR(0), s_rhoL(0), s_rhoR(0), s_pL(1), s_pR(0);
	vector<double> u0(N, uR);
	vector<double> rho0(N, rhoR);
	vector<double> p0(N, pR);
	vector<double> s_u0(N, s_uR);
	vector<double> s_rho0(N,s_rhoR);
	vector<double> s_p0(N,s_pR);
	for (int k=0; k < N*(x_c-xa)/(xb-xa); ++k)
	{
		u0[k] = uL;
		rho0[k] = rhoL;
		p0[k] = pL;
		s_rho0[k] = s_rhoL;
		s_u0[k] = s_uL;
		s_p0[k] = s_pL;
	}
*/	/***************************************/

	/*** Initial and Boundary Conditions ***/
	vector<bool> bc_L(3, false), bc_R(3, false);
	vector<double> VL(6), VR(6); // V = [H, p_tot, p]
	/*** Isentropic transonic ***/
	double rho_init(1.28125), u_init(1.082003561600919), p_init(1.25);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = false;
	/****************************/
	/*** Transonic with shock ***/
/*	double rho_init(1.5), u_init(0.730296743340221), p_init(1.6);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = true;
*/	/****************************/
	/********* Subsonic *********/
/*	double rho_init(1.66875), u_init(0.394721729127866), p_init(1.87);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = true;
*/	/*ifstream if_pstar ("sub_p.dat"); vector<double> pstar(N);
	int k(0);
	while(if_pstar >> pstar[k])
		++k;*/
	/****************************/
	double s_rho_init(0), s_u_init(0), s_p_init(0);
	double gamma(1.4);
	double H_L(4.0), p_tot_L(2), p_L(0);
	double H_R(0), p_tot_R(0), p_R(p_init);
	double s_H_L(0), s_p_tot_L(0), s_p_L(0);
	double s_H_R(0), s_p_tot_R(0), s_p_R(0);
	vector<double> u0(N, u_init);
	vector<double> rho0(N, rho_init);
	vector<double> p0(N, p_init);
	vector<double> s_u0(N, s_u_init);
	vector<double> s_rho0(N,s_rho_init);
	vector<double> s_p0(N,s_p_init);
	VL[0] = H_L;
	VL[1] = p_tot_L;
	VL[2] = p_L;
	VL[3] = s_H_L;
	VL[4] = s_p_tot_L;
	VL[5] = s_p_L;
	VR[0] = H_R;
	VR[1] = p_tot_R;
	VR[2] = p_R;
	VR[3] = s_H_R;
	VR[4] = s_p_tot_R;
	VR[5] = s_p_R;
	/***************************************/

	/******** h and s_h definition *********/
	vector<double> h(N, 1.), dh(N,0.);
	vector<double> h_L(N, 0.), dh_L(N,0.);
	vector<double> h_A(N, 0.), dh_A(N,0.);
	double pi(4*atan(1)), x, xc(0.5), L(0.5), A(1.);
	for (int i = 0; i < N; ++i)
	{
		x = 0.5*dx + i*dx;
		h[i] = 2. - A*(sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
		dh[i] = A*(sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - A*(sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		h_L[i] = - pi/L/L*(x-xc)*A*sin(2*pi*(x-xc)/L)*(x>xc-L/2.)*(x<xc+L/2.);
		dh_L[i] = pi/L/L*(dx*i-xc)*A*sin(2*pi*(dx*i-xc)/L)*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - pi/L/L*(dx*(i+1)-xc)*A*sin(2*pi*(dx*(i+1)-xc)/L)*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		h_A[i] = - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
		dh_A[i] = (sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
	}
	/***************************************/
	
	/************* OPTIMIZATION ************/
	int NP = 2;
	double J(0), J_old(1), coeff(1);
	vector<double> gradJ(NP), param(NP), param_old(NP);
	bool stopcrit (false);
	param[0] = A;
	param[1] = L;
	param_old = param;
	roe_I st_A(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh, h_A, dh_A);
	roe_I st_L(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh, h_L, dh_L);
	vector<roe_I> st;
	st.push_back(st_A);
	st.push_back(st_L);
	vector<vector<vector<double> > > W(NP);
	int time_order (1);
	bool CD (false);
	time_solver TS(t, T, time_order, M, cfl);
	for (int k = 0; k < NP; ++k)
	{
		st[k].set_bc_L(VL, bc_L);
		st[k].set_bc_R(VR, bc_R);
		st[k].set_CD(CD);
		st[k].set_sens_hllc(false);
	}
	unsigned int iter1(0), max_iter1(1000), iter2(0), max_iter2(20);
	double toll = 1e-5;
	TS.solve(st[0], true);
	st[0].get_W(W[0]);
	vector<double> pstar = W[0][2];
	
	param[0] = 0.5;
	param[1] = 0.2;
	A = 0.5; L = 0.2;
	for (int i = 0; i < N; ++i)
	{
		x = 0.5*dx + i*dx;
		h[i] = 2. - A*(sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
		dh[i] = A*(sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - A*(sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		h_L[i] = - pi/L/L*(x-xc)*A*sin(2*pi*(x-xc)/L)*(x>xc-L/2.)*(x<xc+L/2.);
		dh_L[i] = pi/L/L*(dx*i-xc)*A*sin(2*pi*(dx*i-xc)/L)*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - pi/L/L*(dx*(i+1)-xc)*A*sin(2*pi*(dx*(i+1)-xc)/L)*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		h_A[i] = - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
		dh_A[i] = (sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
	}
	for (int k = 0; k < NP; ++k)
	{
		st[k].set_h(h);
		st[k].set_dh(dh);
	}
	st[0].set_sh(h_A);
	st[0].set_dsh(dh_A);
	st[1].set_sh(h_L);
	st[1].set_dsh(dh_L);
	#pragma omp parallel for
	for(int k = 0; k < NP; ++k)
	{
		TS.solve(st[k]);
		st[k].get_W(W[k]);
	}
	J_old = 0.5*L2dot(W[0][2] - pstar, W[0][2] - pstar, dx);
	for (int k = 0; k < NP; ++k)
		fout_param << param[k] << "\t";
	fout_param << endl;
	fout_J << J_old << endl;
	cout << "J = " << J_old << endl;
	while(!stopcrit && iter1 < max_iter1)
	{
		cerr << "Iter " << iter1 << endl;
		++iter1;
		coeff = 1;
		for(int k = 0; k < NP; ++k)
		{
			gradJ[k] = L2dot(W[k][2] - pstar, W[k][5], dx);
			param[k] = param_old[k] - coeff*gradJ[k];
		}
		cerr << "grad_() = ( " << gradJ[0] << ", " << gradJ[1] << ")\n";
		cerr << "New param before checking: " << param[0] << " " << param[1] << endl;
		fout_gradJ << gradJ[0] << "\t" << gradJ[1] << endl;
		/**************** Projection (x_c,L) ****************/
/*		if(!( param[1] > 0 && param[1] < 1 && param[0] < 2*param[1] && param[0] < 2-2*param[1]))
		{
			if (param[0] < 0 && param[1] > 0 && param[1] < 1)
			{
				param[0] = 0.01; // in teoria 0, da cambiare
			}
			if(param[1] > 1 && param[0] < 0.5*param[1]- 0.5)
			{
				param[0] = 0.01;
				param[1] = 1;
			}
			if(param[0] > -2*param[1]+2 && param[0] > 0.5*param[1] - 0.5 && param[0] < 0.5*param[1] + 0.75)
			{
				param[1] = 4./5. + param[1]/5. - 2*param[0]/5.;
				param[0] = -2*param[1]+2;
			}
			if (param[0] > -0.5*param[1] + 5./4. && param[0] > 0.5*param[1]+3./4.)
			{
				param[0] = 1;
				param[1] = 0.5;
			}
			if(param[0] > -0.5*param[1] && param[0] > 2*param[1] && param[0] < -0.5*param[1] + 5./4.)
			{
				param[1] = 2*param[0]/5. + param[1]/5.;
				param[0] = 2*param[1];
			}
			if(param[1] < 0 && param[0] < -0.5*param[1])
			{
				param[0] = 0.01;
				param[1] = 0.01;
			}
			for (int k = 0; k < NP; ++k)
			{
				gradJ[k] = (param_old[k]-param[k])/coeff;
			}
		}
*/		/**********************************************/
		/**************** Projection (A,L) ****************/
		if(!( param[1] > 0 && param[1] < 1 && param[0] < 2 && param[0] > 0))
		{
			if (param[0] < 0)
				param[0] = 0;
			if(param[0] > 2)
				param[0] = 1.99;
			if (param[1] < 0)
				param[1] = 0.01;
			if(param[1]>1)
				param[1] = 1;
			for (int k = 0; k < NP; ++k)
				gradJ[k] = (param_old[k]-param[k])/coeff;
		}
		 /**********************************************/

		A = param[0];
		L = param[1];
		cerr << "New param after checking the domain: " << param[0] << " " << param[1] << endl;
		for (int i = 0; i < N; ++i)
		{
			x = 0.5*dx + i*dx;
			h[i] = 2. - A*(sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
			dh[i] = A*(sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - A*(sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
			h_L[i] = - pi/L/L*(x-xc)*A*sin(2*pi*(x-xc)/L)*(x>xc-L/2.)*(x<xc+L/2.);
			dh_L[i] = pi/L/L*(dx*i-xc)*A*sin(2*pi*(dx*i-xc)/L)*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - pi/L/L*(dx*(i+1)-xc)*A*sin(2*pi*(dx*(i+1)-xc)/L)*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
			h_A[i] = - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
			dh_A[i] = (sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		}
		for (int k = 0; k < NP; ++k)
		{
			st[k].set_h(h);
			st[k].set_dh(dh);
		}
		st[0].set_sh(h_A);
		st[0].set_dsh(dh_A);
		st[1].set_sh(h_L);
		st[1].set_dsh(dh_L);
		
		TS.solve(st[0], true);
		st[0].get_W(W[0]);
		J = 0.5*L2dot(W[0][2] - pstar, W[0][2] - pstar, dx);

		iter2 = 0;
		while(J >= J_old && iter2 < max_iter2)
		{
			if(iter2==0) cerr << "Jold = " << J_old << endl;
			if(iter2 == 0) cerr << "J = " << J << " has not decreased, coeff is halved for the " << iter2+1 << "st time. ";
			if(iter2==1) cerr << "J = " << J << " has not decreased, coeff is halved for the " << iter2+1 << "nd time. ";
			if(iter2==2) cerr << "J = " << J << " has not decreased, coeff is halved for the " << iter2+1 << "rd time. ";
			if(iter2 > 2) cerr << "J = " << J << " has not decreased, coeff is halved for the " << iter2+1 << "th time. ";
			++iter2;
			coeff *= 0.5;
			for (int k = 0; k < NP; ++k)
				param[k] = param_old[k] - coeff*gradJ[k];
			cerr << "Trying with L = " << param[1] << ", A = " << param[0] << endl;
			A = param[0];
			L = param[1];
			for (int i = 0; i < N; ++i)
			{
				x = 0.5*dx + i*dx;
				h[i] = 2. - A*(sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
				dh[i] = A*(sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - A*(sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
			}
			st[0].set_h(h);
			st[0].set_dh(dh);
			TS.solve(st[0], true);
			st[0].get_W(W[0]);
			J = 0.5*L2dot(W[0][2]-pstar, W[0][2]-pstar, dx);
		}
		if(iter2 == max_iter2)
			cerr << "Max_iter2 raggiunto.\n";
		stopcrit = fabs(param_old[0]-param[0]) < toll;
		for (int k = 1; k < NP; ++k)
			stopcrit = stopcrit && (fabs(param_old[k]-param[k]) < toll);
		cerr << "After second while, stopcrit = " << stopcrit << endl;
		if(J < J_old)
		{
			for (int i = 0; i < N; ++i)
			{
				x = 0.5*dx + i*dx;
				h_L[i] = - pi/L/L*(x-xc)*A*sin(2*pi*(x-xc)/L)*(x>xc-L/2.)*(x<xc+L/2.);
				dh_L[i] = pi/L/L*(dx*i-xc)*A*sin(2*pi*(dx*i-xc)/L)*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - pi/L/L*(dx*(i+1)-xc)*A*sin(2*pi*(dx*(i+1)-xc)/L)*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
				h_A[i] = - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
				dh_A[i] = (sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
			}

			st[0].set_sh(h_A);
			st[0].set_dsh(dh_A);
			st[1].set_sh(h_L);
			st[1].set_dsh(dh_L);
			st[1].set_h(h);
			st[1].set_dh(dh);			
			#pragma omp parallel for
			for(int k = 0; k < NP; ++k)
			{
				TS.solve(st[k]);
				st[k].get_W(W[k]);
			}
		}
		
		param_old = param;
		J_old = J;
		for (int k = 0; k < NP; ++k)
			fout_param << param[k] << "\t";
		fout_param << endl;
		fout_J << J << endl;
		cout << "With param = (" << A << ", " << L << ")\t J = " << J << endl;
	}
	for(int k = 0; k < NP; ++k)
	{
		st[k].print_physical(path+to_string(k)+"_");
	}
	/***************************************/
	return 0;
}
