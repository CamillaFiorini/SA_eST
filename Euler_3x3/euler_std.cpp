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
	double xa(0), xb(1), dx(1e-3), T(0.1), t(0), cfl(0.5);
	mesh M (xa, xb, dx);
	string path = "results/"; //"../../results/Euler_3x3_q1d/err_extrapol/isentropic/diff_ord1/dx5e-3/big_da/da005/";
	double da(0);
	ofstream file_da(path+"da.dat");
	file_da << da << endl;
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
	vector<vector<double> > W;
	vector<double> h(N, 1.), dh(N,0.);
	roe_I st(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh);
	st.set_bc_L(VL, bc_L);
	st.set_bc_R(VR, bc_R);
	int time_order (1);
	bool CD (false);
	st.set_CD(CD);
	st.set_sens_hllc(false);
	st.print_physical(path);
	time_solver TS(t, T, time_order, M, cfl);
	TS.solve(st);
	st.print_physical(path, ios::out | ios::app);
	return 0;
}
