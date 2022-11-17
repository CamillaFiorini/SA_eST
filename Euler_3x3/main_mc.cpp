#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include<chrono>
#include<sstream>
#include<string>
#include<random>
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
	double xa(0), xb(1), dx(1e-3), T(0.2), t(0), cfl(0.5);
	mesh M (xa, xb, dx);
	string path = "new_sub_Euler1d_CD/";
	int N = M.get_N();
	/***************************************/
	/************* Shock tube **************/
	double uL(0), uR(0), rhoL(1), rhoR(0.125), pL(1.0), pR(0.1), x_c(0.5), gamma(1.4);
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
	/***************************************/
	/*** Initial and Boundary Conditions ***/
	vector<bool> bc_L(3, false), bc_R(3, false);
	vector<double> VL(6), VR(6); // V = [H, p_tot, p]
	/*** Isentropic transonic ***/
/*	double rho_init(1.28125), u_init(1.082003561600919), p_init(1.25);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = false;
*/	/****************************/
	/*** Transonic with shock ***/
/*	double rho_init(1.5), u_init(0.730296743340221), p_init(1.6);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = true;
*/	/****************************/
	/********* Subsonic *********/
/*	double rho_init(1.66875), u_init(0.394721729127866), p_init(1.87);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = true;
*/	/****************************/
/*	double s_rho_init(0), s_u_init(0), s_p_init(0);
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
*/	/***************************************/

	/******** h and s_h definition *********/
	vector<double> h(N, 1.), dh(N,0.);
//	double pi(4*atan(1)), x, xc(0.5), L(0.5);
	/***************************************/
	
	/************** UQ intro ***************/
//	roe_I st(rho0,u0, p0, s_rho0,s_u0, s_p0, gamma, h, dh);
	vector<vector<double> > W(6), average(3), variance(3), param(6);
	int time_order (1);
	bool CD (true);
	time_solver TS(t, T, time_order, M, cfl);
/*	st.set_bc_L(VL, bc_L);
	st.set_bc_R(VR, bc_R);
	st.set_CD(CD);
	st.set_sens_hllc(false);
*/	int MC_max = 1000;
	/***************************************/
	ofstream rho_ave_out (path+"rho_ave.dat");
	ofstream u_ave_out (path+"u_ave.dat");
	ofstream p_ave_out (path+"p_ave.dat");
	ofstream rho_var_out (path+"rho_var.dat");
	ofstream u_var_out (path+"u_var.dat");
	ofstream p_var_out (path+"p_var.dat");
	ofstream param_out (path+"param.dat");
	ofstream hist_rho_out (path+"hist_rho.dat");
	ofstream hist_u_out (path+"hist_u.dat");
	ofstream hist_p_out (path+"hist_p.dat");
	/***************************************/
	
	default_random_engine generator;
/*	normal_distribution<double> distribution_L(L,sqrt(L*0.001));
	normal_distribution<double> distribution_xc(xc,sqrt(xc*0.001));
	normal_distribution<double> distribution_H_L(H_L,sqrt(H_L*0.001));
	normal_distribution<double> distribution_p_tot_L(p_tot_L,sqrt(p_tot_L*0.001));
	normal_distribution<double> distribution_p_R(p_R,sqrt(p_R*0.001));
*/
	normal_distribution<double> distribution_rhoL(rhoL,sqrt(rhoL*0.001));
	normal_distribution<double> distribution_rhoR(rhoR,sqrt(rhoR*0.001));
	normal_distribution<double> distribution_uL(uL,sqrt(pR*0.001));
	normal_distribution<double> distribution_uR(uR,sqrt(pR*0.001));
	normal_distribution<double> distribution_pL(pL,sqrt(pL*0.001));
	normal_distribution<double> distribution_pR(pR,sqrt(pR*0.001));
	
	for (int i = 0; i < 3; ++i)
	{
		average[i].assign(N, 0);
		variance[i].assign(N, 0);
	}
	for (int i = 0; i < 6; ++i)
		param[i].resize(MC_max);
/*	ifstream if_param ("param.dat");
	double dummy; int k = 0;
	string str;
	while(getline(if_param, str))
	{
		stringstream iss (str);
		int i = 0;
		while(iss >> dummy)
		{	
			param[i][k] = dummy;
			++i;
		}
		++k;
	}	
*/
	for(int iter_MC = 0; iter_MC < MC_max; ++iter_MC)
	{
		cout << iter_MC << endl;
/*		L = distribution_L(generator);
		xc = distribution_xc(generator);
		H_L = distribution_H_L(generator);
		p_tot_L = distribution_p_tot_L(generator);
		p_R = distribution_p_R(generator);
		param[0][iter_MC] = L;
		param[1][iter_MC] = xc;
		param[2][iter_MC] = H_L;
		param[3][iter_MC] = p_tot_L;
		param[4][iter_MC] = p_R;
 */
		rhoL = distribution_rhoL(generator);
		rhoR = distribution_rhoR(generator);
		uL = distribution_uL(generator);
		uR = distribution_uR(generator);
		pL = distribution_pL(generator);
		pR = distribution_pR(generator);
		param[0][iter_MC] = rhoL;
		param[1][iter_MC] = rhoR;
		param[2][iter_MC] = uL;
		param[3][iter_MC] = uR;
		param[4][iter_MC] = pL;
		param[5][iter_MC] = pR;
		
		u0.assign(N, uR);
		rho0.assign(N, rhoR);
		p0.assign(N, pR);
		s_u0.assign(N, s_uR);
		s_rho0.assign(N,s_rhoR);
		s_p0.assign(N,s_pR);
		for (int k=0; k < N*(x_c-xa)/(xb-xa); ++k)
		{
			u0[k] = uL;
			rho0[k] = rhoL;
			p0[k] = pL;
			s_rho0[k] = s_rhoL;
			s_u0[k] = s_uL;
			s_p0[k] = s_pL;
		}
		
		for (int i = 0; i < 6; ++i)
			param_out << param[i][iter_MC] << "\t";
		param_out << endl;
		
/*		for (int i = 0; i < N; ++i)
		{
			x = 0.5*dx + i*dx;
			h[i] = 2. - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
			dh[i] = (sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		}
		st.set_h(h);
		st.set_dh(dh);
		VL[0] = H_L;
		VL[1] = p_tot_L;
		VR[2] = p_R;
		st.set_bc_L(VL, bc_L);
		st.set_bc_R(VR, bc_R);
*/
		roe_I st(rho0,u0, p0, s_rho0,s_u0, s_p0, gamma, h, dh);
		st.set_sens_hllc(false);
		st.set_CD(CD);
		TS.solve(st, true);
		st.get_W(W);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < N; ++j)
				average[i][j] += W[i][j];
		if (iter_MC < 10)
		{
			st.print_physical(path, ios::out | ios::app);
		}
		// addition for resumbission
		hist_rho_out << W[0][350] << "\t" << W[0][600] << "\t" << W[0][850] << endl;
		hist_u_out << W[1][350] << "\t" << W[1][600] << "\t" << W[1][850] << endl;
		hist_p_out << W[2][350] << "\t" << W[2][600] << "\t" << W[2][850] << endl;
	}

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < N; ++j)
			average[i][j] /= MC_max;
	cerr << "Printing average on file: ";
	for (int j = 0; j < N; ++j)
	{
		//cerr << average[0][j] << "\t" << average[1][j] << "\t" << average[2][j] << endl;
		rho_ave_out << average[0][j] << "\t";
		u_ave_out << average[1][j] << "\t";
		p_ave_out << average[2][j] << "\t";
	}
	cerr << "\nAverage printed" << endl;

	for(int iter_MC = 0; iter_MC < MC_max; ++iter_MC)
	{
		cout << iter_MC << endl;
/*		L = param[0][iter_MC];
		xc = param[1][iter_MC];
		H_L = param[2][iter_MC];
		p_tot_L = param[3][iter_MC];
		p_R = param[4][iter_MC];
		for (int i = 0; i < N; ++i)
		{
			x = 0.5*dx + i*dx;
			h[i] = 2. - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
			dh[i] = (sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		}
		st.set_h(h);
		st.set_dh(dh);
		VL[0] = H_L;
		VL[1] = p_tot_L;
		VR[2] = p_R;
		st.set_bc_L(VL, bc_L);
		st.set_bc_R(VR, bc_R);
	*/
		rhoL = param[0][iter_MC];
		rhoR = param[1][iter_MC];
		uL = param[2][iter_MC];
		uR = param[3][iter_MC];
		pL = param[4][iter_MC];
		pR = param[5][iter_MC];
		
		u0.assign(N, uR);
		rho0.assign(N, rhoR);
		p0.assign(N, pR);
		s_u0.assign(N, s_uR);
		s_rho0.assign(N,s_rhoR);
		s_p0.assign(N,s_pR);
		for (int k=0; k < N*(x_c-xa)/(xb-xa); ++k)
		{
			u0[k] = uL;
			rho0[k] = rhoL;
			p0[k] = pL;
			s_rho0[k] = s_rhoL;
			s_u0[k] = s_uL;
			s_p0[k] = s_pL;
		}
		roe_I st(rho0,u0, p0, s_rho0,s_u0, s_p0, gamma, h, dh);
		st.set_sens_hllc(false);
		st.set_CD(CD);
		TS.solve(st, true);
		st.get_W(W);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < N; ++j)
				variance[i][j] += (W[i][j] - average[i][j])*(W[i][j] - average[i][j]);
	}
	
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < N; ++j)
			variance[i][j] /= MC_max-1;
	for (int j = 0; j < N; ++j)
	{
		rho_var_out << variance[0][j] << "\t";
		u_var_out << variance[1][j] << "\t";
		p_var_out << variance[2][j] << "\t";
	}
	/***************************************/

	return 0;
}
