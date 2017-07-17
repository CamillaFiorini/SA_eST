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
	double xa(0), xb(1), dx(1e-3), T(0.1), t(0);
	//int n_print(5/(100*dx));
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(0), uR(0), rhoL(1), rhoR(0.125), pL(1), pR(0.1), x_c(0.5), gamma(1.4);//uL(19.5975), uR(-6.19633), rhoL(5.99924), rhoR(5.99242), pR(46.0950), pL(460.894), x_c(0.5), gamma(1.4);
	double s_uL(0), s_uR(0), s_rhoL(0), s_rhoR(0), s_pL(1), s_pR(0);
	double cfl(0.5);
	vector<double> u0(N, uR);
	vector<double> rho0(N, rhoR);
	vector<double> p0(N, pR);
	vector<double> s_u0(N, s_uR);
	vector<double> s_rho0(N,s_rhoR);
	vector<double> s_p0(N,s_pR);
	cout.precision(15);
	for (int k=0; k < N*x_c; ++k)
	{
		u0[k] = uL;
		rho0[k] = rhoL;
		p0[k] = pL;
		s_rho0[k] = s_rhoL;
		s_u0[k] = s_uL;
		s_p0[k] = s_pL;
	}
	
	vector<vector<double> > W;
	
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
	roe_II st(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma);
	int time_order (2);
	bool CD (true);
	st.set_CD(CD);
	st.set_sens_hllc(false);
	time_solver TS(t, T, time_order, M, cfl);
	/*
	ofstream file_u ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/u.dat");
	ofstream file_rho ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/rho.dat");
	ofstream file_p ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/p.dat");
	ofstream file_s_u ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/s_u.dat");
	ofstream file_s_rho ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/s_rho.dat");
	ofstream file_s_p ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/s_p.dat");
	ofstream file_t ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/t.dat");
	ofstream file_d1 ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/d1.dat");
	ofstream file_d2 ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/d2.dat");
	ofstream file_c ("../../results/Euler_3x3_conv/second_order/Diff_HLL/dx1e-5/c.dat");
	*/
	ofstream file_u ("results/u.dat");
	ofstream file_rho ("results/rho.dat");
	ofstream file_p ("results/p.dat");
	ofstream file_s_u ("results/s_u.dat");
	ofstream file_s_rho ("results/s_rho.dat");
	ofstream file_s_p ("results/s_p.dat");
	ofstream file_t ("results/t.dat");
	ofstream file_d1 ("results/d1.dat");
	ofstream file_d2 ("results/d2.dat");
	ofstream file_c ("results/c.dat");
	
	file_u.precision(15);
	file_rho.precision(15);
	file_p.precision(15);
	file_s_u.precision(15);
	file_s_rho.precision(15);
	file_s_p.precision(15);
	file_t.precision(15);
	
	for (int k = 0; k < N; ++k)
	{
		file_u << u0[k] << "\t";
		file_rho << rho0[k] << "\t";
		file_p << p0[k] << "\t";
		file_s_u << s_u0[k] << "\t";
		file_s_rho << s_rho0[k] << "\t";
		file_s_p << s_p0[k] << "\t";
	}
	file_u << endl;
	file_rho << endl;
	file_p << endl;
	file_s_rho << endl;
	file_s_u << endl;
	file_s_p << endl;
	file_t << t << endl;
	
	TS.solve(st);
	
	st.get_W(W);
	for (int k=0; k<N; ++k)
	{
		file_rho << W[0][k] << "\t";
		file_u << W[1][k] << "\t";
		file_p << W[2][k] << "\t";
		file_s_rho << W[3][k] << "\t";
		file_s_u << W[4][k] << "\t";
		file_s_p << W[5][k] << "\t";
	}
	file_rho << endl;
	file_u << endl;
	file_p << endl;
	file_s_rho << endl;
	file_s_u << endl;
	file_s_p << endl;
	file_t << T << endl;

	return 0;
}
