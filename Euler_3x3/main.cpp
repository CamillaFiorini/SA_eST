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
	double xa(0), xb(1), dx(5e-3), T(0.1), t(0);
	//int n_print(5/(100*dx));
	string path = "results/";
	//string path = "../../results/Euler_3x3_conv/second_order/Diff_HLLC/dx1e-2/";
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(0), uR(0), rhoL(1), rhoR(0.125), pL(1), pR(0.1), x_c(0.5), gamma(1.4);
	double s_uL(0), s_uR(0), s_rhoL(0), s_rhoR(0), s_pL(1), s_pR(0);
	double cfl(0.5);
	vector<double> u0(N, uR);
	vector<double> rho0(N, rhoR);
	vector<double> p0(N, pR);
	vector<double> s_u0(N, s_uR);
	vector<double> s_rho0(N,s_rhoR);
	vector<double> s_p0(N,s_pR);
	vector<double> h(N, 1.), dh(N,0.);
	double pi(4*atan(1)), x, L(0.25);
	for (int i = 0; i < N; ++i)
	{
		x = 0.5*dx + i*dx;
		h[i] = 2. - (sin((x-x_c)/L*pi - 0.5*pi)*sin((x-x_c)/L*pi - 0.5*pi))*(x>x_c-L/2.)*(x<x_c+L/2.);//2-exp(-(x-x_c)*(x-x_c)/0.004);
		dh[i] = //-2*dx*sin((x-x_c)/L*pi - 0.5*pi)*cos((x-x_c)/L*pi - 0.5*pi)*pi/L*(x>x_c-L/2.)*(x<x_c+L/2.);
		(sin((dx*i-x_c)/L*pi - 0.5*pi)*sin((dx*i-x_c)/L*pi - 0.5*pi))*(dx*i>x_c-L/2.)*(dx*i<x_c+L/2.) - (sin((dx*(i+1)-x_c)/L*pi - 0.5*pi)*sin((dx*(i+1)-x_c)/L*pi - 0.5*pi))*(dx*(i+1)>x_c-L/2.)*(dx*(i+1)<x_c+L/2.);
	}
	cout.precision(15);
	
	for (int k=0; k < N*(x_c-xa)/(xb-xa); ++k)
	{
		u0[k] = uL;
		rho0[k] = rhoL;
		p0[k] = pL;
		s_rho0[k] = s_rhoL;
		s_u0[k] = s_uL;
		s_p0[k] = s_pL;
	}
	/********** Restoring ***************/
/*	t = 0.0114115027477935;
	ifstream if_rho ("rho_int.txt");
	ifstream if_p ("p_int.txt");
	ifstream if_u ("u_int.txt");
	ifstream if_s_rho ("s_rho_int.txt");
	ifstream if_s_p ("s_p_int.txt");
	ifstream if_s_u ("s_u_int.txt");
	double dummy; int k(0);
	while(if_rho >> dummy)
	{
		rho0[k] = dummy;
		++k;
	}
	k = 0;
	while(if_p >> dummy)
	{
		p0[k] = dummy;
		++k;
	}
	k = 0;
	while(if_u >> dummy)
	{
		u0[k] = dummy;
		++k;
	}
	k = 0;
	while(if_s_rho >> dummy)
	{
		s_rho0[k] = dummy;
		++k;
	}
	k = 0;
	while(if_s_p >> dummy)
	{
		s_p0[k] = dummy;
		++k;
	}
	k = 0;
	while(if_s_u >> dummy)
	{
		s_u0[k] = dummy;
		++k;
	}
	cout << "End of restoring\n"; */
	/****************************************/
	roe_I st(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma, h, dh);
	int time_order (1);
	bool CD (false);
	st.set_CD(CD);
	st.set_sens_hllc(false);
	time_solver TS(t, T, time_order, M, cfl);
	
	st.print_physical(path);
	cout << TS.solve(st) << endl;
	st.print_physical(path, ios::out | ios::app);
	
	return 0;
}
