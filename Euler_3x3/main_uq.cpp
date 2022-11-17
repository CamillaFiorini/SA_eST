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
	double xa(0), xb(1), dx(1e-3), T(0.2), t(0), cfl(0.5);
	mesh M (xa, xb, dx);
	string path = "res_uq/RiemannPb/"; //"../../results/Euler_3x3_q1d/err_extrapol/isentropic/diff_ord1/dx5e-3/big_da/da005/";
	int N = M.get_N();
	/***************************************/

	/*** Initial and Boundary Conditions ***/
	int NP = 6;
	vector<bool> bc_L(3, false), bc_R(3, false);
	vector<vector<double> > VL(NP), VR(NP); // V = [H, p_tot, p]
	for (int k = 0; k < NP; ++k)
	{
		VL[k].resize(6);
		VR[k].resize(6);
	}
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
	double H_L(4), p_tot_L(2), p_L(0);
	double H_R(0), p_tot_R(0), p_R(p_init);
	double s_H_L(0), s_p_tot_L(0), s_p_L(0);
	double s_H_R(0), s_p_tot_R(0), s_p_R(0);
	vector<double> u0(N, u_init);
	vector<double> rho0(N, rho_init);
	vector<double> p0(N, p_init);
	vector<double> s_u0(N, s_u_init);
	vector<double> s_rho0(N,s_rho_init);
	vector<double> s_p0(N,s_p_init);
	for(int k = 0; k < NP; ++k)
	{
		VL[k][0] = H_L;
		VL[k][1] = p_tot_L;
		VL[k][2] = p_L;
		VL[k][3] = s_H_L;
		VL[k][4] = s_p_tot_L;
		VL[k][5] = s_p_L;
		VR[k][0] = H_R;
		VR[k][1] = p_tot_R;
		VR[k][2] = p_R;
		VR[k][3] = s_H_R;
		VR[k][4] = s_p_tot_R;
		VR[k][5] = s_p_R;
	}
	VL[2][3] = 1; // s_H_L = 1
	VL[3][4] = 1; // s_p_tot_L = 1
	VR[4][5] = 1; // s_p_R = 1
*/	/***************************************/

	/******** h and s_h definition *********/
	vector<double> h(N, 1.), dh(N,0.);
/*	vector<double> h_xc(N, 0.), dh_xc(N,0.);
	vector<double> h_L(N, 0.), dh_L(N,0.);
	double pi(4*atan(1)), x, xc(0.5), L(0.5);
	for (int i = 0; i < N; ++i)
	{
		x = 0.5*dx + i*dx;
		h[i] = 2. - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);
		dh[i] = (sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		h_xc[i] = - pi/L*sin(2*pi*(x-xc)/L)*(x>xc-L/2.)*(x<xc+L/2.);
		dh_xc[i] = pi/L*sin(2*pi*(dx*i-xc)/L)*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - pi/L*sin(2*pi*(dx*(i+1)-xc)/L)*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
		h_L[i] = - pi/L/L*(x-xc)*sin(2*pi*(x-xc)/L)*(x>xc-L/2.)*(x<xc+L/2.);
		dh_L[i] = pi/L/L*(dx*i-xc)*sin(2*pi*(dx*i-xc)/L)*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - pi/L/L*(dx*(i+1)-xc)*sin(2*pi*(dx*(i+1)-xc)/L)*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
	}
*/	/***************************************/
	
	/******* For the Riemann problem *******/
	double uL(0), uR(0), rhoL(1), rhoR(0.125), pL(1.0), pR(0.1), gamma(1.4);
	vector<double> u0(N, uR);
	vector<double> rho0(N, rhoR);
	vector<double> p0(N, pR);
	for (int i = 0; i < N/2; ++i)
	{
		rho0[i] = rhoL;
		u0[i] = uL;
		p0[i] = pL;
	}
	vector<vector<double> > s_u0(NP, u0);
	vector<vector<double> > s_rho0(NP,u0);
	vector<vector<double> > s_p0(NP,u0);
	for (int i = 0; i < N/2; ++i)
	{
		s_rho0[0][i] = 1;
		s_rho0[1][i+N/2] = 1;
		s_u0[2][i] = 1;
		s_u0[3][i+N/2] = 1;
		s_p0[4][i] = 1;
		s_p0[5][i+N/2] = 1;
	}
	int time_order (1);
	bool CD (true);
	time_solver TS(t, T, time_order, M, cfl);
	vector<roe_I> st;
	for (int k = 0; k < NP; ++k)
	{
		roe_I dummy(rho0,u0, p0, s_rho0[k],s_u0[k], s_p0[k],gamma, h, dh);
		st.push_back(dummy);
		st[k].set_bc_L(VL[k], bc_L);
		st[k].set_bc_R(VR[k], bc_R);
		st[k].set_CD(CD);
		st[k].set_sens_hllc(false);
	}
	/***************************************/
	
	/****************** UQ *****************/
/*	roe_I st_L(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh, h_L, dh_L);
	roe_I st_xc(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh, h_xc, dh_xc);
	roe_I st_bc(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh);
	vector<roe_I> st;
	st.push_back(st_L);
	st.push_back(st_xc);
	st.push_back(st_bc);
	st.push_back(st_bc);
	st.push_back(st_bc);
	int time_order (1);
	bool CD (false);
	time_solver TS(t, T, time_order, M, cfl);
	for (int k = 0; k < NP; ++k)
	{
		st[k].set_bc_L(VL[k], bc_L);
		st[k].set_bc_R(VR[k], bc_R);
		st[k].set_CD(CD);
		st[k].set_sens_hllc(false);
	}*/
	
	
	#pragma omp parallel for
	for(int k = 0; k < NP; ++k)
	{
		printf("thread boh, iter %d\n", k);
		st[k].print_physical(path+"CD_a"+to_string(k)+"_");
		TS.solve(st[k]);
		st[k].print_physical(path+"CD_a"+to_string(k)+"_", ios::out | ios::app);
	}
	/***************************************/
	return 0;
}
