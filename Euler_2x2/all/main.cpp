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
	double xa(0), xb(1), dx(1e-3), T(0.1), t(0);
	mesh M (xa, xb, dx);
	string path("results/isolated_shock/");
	int N (M.get_N());
	double uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4), x_c(0.5); // uL(0), uR(0), tauL(0.7), tauR(0.4), gamma(1.4), x_c(0.5);//
	double cfl(0.5);
	vector<double> u0(N, uR);
	vector<double> tau0(N, tauR);
	vector<double> s_u0(N, 0);
	vector<double> s_tau0(N,0);
	cout.precision(15);

	for (int k=0; k < N*x_c; ++k)
	{
		u0[k] = uL;
		tau0[k] = tauL;
		s_tau0[k] = 1; // 0;//
		s_u0[k] = -9.351212140372281;//0;//
	}
	/************* Bump *************
	double pi (4*atan(1)), L(0.25), m(0.05);
	for (int k=0; k < N; ++k)
	{
		double x = (k+0.5)*dx;
		if(x > x_c - 0.5*L && x < x_c)
		{
			tau0[k] = tauL+(m-tauL)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[k] = pi/L*sin(2*pi/L*(x-x_c));//-sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
		}
		if(x > x_c && x < x_c + 0.5*L)
		{
			tau0[k] = tauR+(m-tauR)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			s_tau0[k] = pi/L*sin(2*pi/L*(x-x_c));//-sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
		}
	}
	********************************/
	
	roe_I st(tau0,u0,s_tau0,s_u0,gamma);
	int time_order (1);
	bool CD_state (true), CD_sens(true);
    st.set_CD_state(CD_state);
    st.set_CD_sens(CD_sens);
	st.print_physical(path);
	time_solver TS(t, T, time_order, M, cfl, path);
	TS.solve(st);
	st.print_physical(path, ios::out | ios::app);

	return 0;
}
