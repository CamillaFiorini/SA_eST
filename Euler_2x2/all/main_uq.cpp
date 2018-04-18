#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include<chrono>
#include<string>
#include"mesh.hpp"
#include"state.hpp"
#include"godunov.hpp"
#include"roe.hpp"
#include"roe_I.hpp"
#include"roe_II.hpp"
#include"utilities.hpp"
#include"time_solver.hpp"

using namespace std;

int main()
{
	double xa(0), xb(1), dx(2e-3), T(0.03), t(0);
	mesh M (xa, xb, dx);
	string path("res_uq/");
	int N (M.get_N());
	double uL(0), uR(0), tauL(0.7), tauR(0.4), gamma(1.4), x_c(0.5);//uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4), x_c(0.5);
	double cfl(0.5);
	int NP = 5;
	vector<double> u0(N, uR);
	vector<double> tau0(N, tauR);
	vector<double> s_u0(N, 0);
	vector<vector<double> > s_tau0(NP,u0);
	cout.precision(15);

	for (int k=0; k < N*x_c; ++k)
	{
		u0[k] = uL;
		tau0[k] = tauL;
		s_u0[k] = 0;
	}
	/************* Bump *************/
	double pi (4*atan(1)), L(0.25), m(0.05);
	for (int k=0; k < N; ++k)
	{
		double x = (k+0.5)*dx;
		if(x <= x_c - 0.5*L)
		{
			s_tau0[0][k] = 1;
			s_tau0[1][k] = 0;
			s_tau0[2][k] = 0;
			s_tau0[3][k] = 0;
			s_tau0[4][k] = 0;

		}
		if(x > x_c - 0.5*L && x <= x_c)
		{
			tau0[k] = tauL+(m-tauL)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[0][k] = 1-sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[1][k] = 0;
			s_tau0[2][k] = sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[3][k] = (m-tauL)*x_c*pi/L/L*sin(2*pi/L*(x-x_c)+pi);
			s_tau0[4][k] = (m-tauL)*(-pi)/L*sin(2*pi/L*(x-x_c)+pi);
		}
		if(x > x_c && x <= x_c + 0.5*L)
		{
			tau0[k] = tauR+(m-tauR)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[0][k] = 0;
			s_tau0[1][k] = 1-sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[2][k] = sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[3][k] = (m-tauR)*x_c*pi/L/L*sin(2*pi/L*(x-x_c)+pi);
			s_tau0[4][k] = (m-tauR)*(-pi)/L*sin(2*pi/L*(x-x_c)+pi);
		}
		if(x > x_c + 0.5*L)
		{
			s_tau0[0][k] = 0;
			s_tau0[1][k] = 1;
			s_tau0[2][k] = 0;
			s_tau0[3][k] = 0;
			s_tau0[4][k] = 0;
		}
	}
	/********************************/
	
	roe_II st_tauL(tau0,u0,s_tau0[0],s_u0,gamma);
	roe_II st_tauR(tau0,u0,s_tau0[1],s_u0,gamma);
	roe_II st_m(tau0,u0,s_tau0[2],s_u0,gamma);
	roe_II st_L(tau0,u0,s_tau0[3],s_u0,gamma);
	roe_II st_xc(tau0,u0,s_tau0[4],s_u0,gamma);
	vector<roe_II> st;
	st.push_back(st_tauL);
	st.push_back(st_tauR);
	st.push_back(st_m);
	st.push_back(st_L);
	st.push_back(st_xc);
	int time_order (2);
	bool CD (true);
	for(int k = 0; k < NP; ++k)
		st[k].set_CD(CD);
	time_solver TS(t, T, time_order, M, cfl);
	
	/*st_tauL.print_physical(path+"a0_");
	TS.solve(st_tauL);
	st_tauL.print_physical(path+"a0_", ios::out | ios::app);*/
	
	for(int k = 0; k < NP; ++k)
	{
		cerr << k << "\tst[k].size() = " << st.size() << endl;
		st[k].print_physical(path+"a"+to_string(k)+"_");
		TS.solve(st[k]);
		st[k].print_physical(path+"a"+to_string(k)+"_", ios::out | ios::app);
	}
	return 0;
}