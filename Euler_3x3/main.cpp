#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include<chrono>
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
	double xa(0), xb(1), dx(1e-3), T(100), t(0), cfl(0.5);
	mesh M (xa, xb, dx);
	string path = "results/"; //"../../results/Euler_3x3_q1d/err_extrapol/isentropic/diff_ord1/dx5e-3/big_da/da005/";
	double da(0);
	ofstream file_da(path+"da.dat");
	file_da << da << endl;
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
/*	double rho_init(1.28125), u_init(1.082003561600919), p_init(1.25);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = false;
*/	/****************************/
	/*** Transonic with shock ***/
	double rho_init(1.5), u_init(0.730296743340221), p_init(1.6);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = true;
	/****************************/
	/********* Subsonic *********/
/*	double rho_init(1.66875), u_init(0.394721729127866), p_init(1.87);
	bc_L[0] = true; bc_L[1] = true; bc_L[2] = false;
	bc_R[0] = false; bc_R[1] = false; bc_R[2] = true;
*/	/****************************/
	double s_rho_init(0), s_u_init(0), s_p_init(0);
	double gamma(1.4);
	double H_L(4.0+da), p_tot_L(2), p_L(0);
	double H_R(0), p_tot_R(0), p_R(p_init);
	double s_H_L(1), s_p_tot_L(0), s_p_L(0);
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

	/************ h definition *************/
	vector<double> h(N, 1.), dh(N,0.);
	double pi(4*atan(1)), x, xc(0.5), L(0.5);
	for (int i = 0; i < N; ++i)
	{
		x = 0.5*dx + i*dx;
		h[i] = 2. - (sin((x-xc)/L*pi - 0.5*pi)*sin((x-xc)/L*pi - 0.5*pi))*(x>xc-L/2.)*(x<xc+L/2.);//2-exp(-(x-xc)*(x-xc)/0.004);
		dh[i] = //-2*dx*sin((x-xc)/L*pi - 0.5*pi)*cos((x-xc)/L*pi - 0.5*pi)*pi/L*(x>xc-L/2.)*(x<xc+L/2.);
		(sin((dx*i-xc)/L*pi - 0.5*pi)*sin((dx*i-xc)/L*pi - 0.5*pi))*(dx*i>xc-L/2.)*(dx*i<xc+L/2.) - (sin((dx*(i+1)-xc)/L*pi - 0.5*pi)*sin((dx*(i+1)-xc)/L*pi - 0.5*pi))*(dx*(i+1)>xc-L/2.)*(dx*(i+1)<xc+L/2.);
	}
	/***************************************/
	
	roe_I st(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh);
	st.set_bc_L(VL, bc_L);
	st.set_bc_R(VR, bc_R);
	int time_order (1);
	bool CD (false);
	st.set_CD(CD);
	st.set_sens_hllc(false);
	time_solver TS(t, T, time_order, M, cfl);
	st.print_physical(path);
	cerr << TS.solve(st) << endl;
	st.print_physical(path, ios::out | ios::app);
	
	return 0;
}
