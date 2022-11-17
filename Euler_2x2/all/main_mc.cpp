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
#include"godunov.hpp"
#include"roe.hpp"
#include"roe_I.hpp"
#include"roe_II.hpp"
#include"utilities.hpp"
#include"time_solver.hpp"

using namespace std;

int main()
{
	/********* Domain definition ***********/
	double xa(0), xb(1), dx(2e-3), T(0.03), t(0), cfl(0.5);
	mesh M (xa, xb, dx);
	string path = "res_mc/riemann_pb/";
	int N = M.get_N();
	/***************************************/
	/**************** IC def ***************/
	double tauL(0.7), tauR(0.2), uR(0), uL(0), gamma(1.4), x_c(0.5);
	//double tauL(0.7), tauR(0.4), gamma(1.4), x_c(0.5), pi (4*atan(1)), L(0.25), m(0.05);
	vector<double> u0(N, 0);
	vector<double> tau0(N, 0);
	vector<double> s_u0(N, 0);
	vector<double> s_tau0(N,0);
	/***************************************/
	/************** MC intro ***************/
	int D (4), NP(4);
	vector<vector<double> > W(D), average(D/2), variance(D/2), param(NP);
	int time_order (2);
	bool CD (true);
	time_solver TS(t, T, time_order, M, cfl);
	int MC_max = 500;
	/***************************************/
	ofstream tau_ave_out (path+"tau_ave.dat");
	ofstream u_ave_out (path+"u_ave.dat");
	ofstream tau_var_out (path+"tau_var.dat");
	ofstream u_var_out (path+"u_var.dat");
	ofstream param_out (path+"param.dat");
	/***************************************/
	
	default_random_engine generator;
/*
	normal_distribution<double> distribution_tauL(tauL,sqrt(tauL*0.001));
	normal_distribution<double> distribution_tauR(tauR,sqrt(tauR*0.001));
	normal_distribution<double> distribution_m(m,sqrt(m*0.001));
	normal_distribution<double> distribution_L(L,sqrt(L*0.001));
	normal_distribution<double> distribution_xc(x_c,sqrt(x_c*0.001));
*/
	normal_distribution<double> distribution_tauL(tauL,sqrt(tauL*0.001));
	normal_distribution<double> distribution_tauR(tauR,sqrt(tauR*0.001));
	normal_distribution<double> distribution_uL(uL,sqrt(0.0001));
	normal_distribution<double> distribution_uR(uR,sqrt(0.0001));

	for (int i = 0; i < D/2; ++i)
	{
		average[i].assign(N, 0);
		variance[i].assign(N, 0);
	}
	for (int i = 0; i < NP; ++i)
		param[i].resize(MC_max);

	for(int iter_MC = 0; iter_MC < MC_max; ++iter_MC)
	{
		cout << iter_MC << endl;
/*
		tauL = distribution_tauL(generator);
		tauR = distribution_tauR(generator);
		m = distribution_m(generator);
		L = distribution_L(generator);
		x_c = distribution_xc(generator);
		param[0][iter_MC] = tauL;
		param[1][iter_MC] = tauR;
		param[2][iter_MC] = m;
		param[3][iter_MC] = L;
		param[4][iter_MC] = x_c;
		for (int k = 0; k < N; ++k)
		{
			double x = (k+0.5)*dx;
			if(x <= x_c - 0.5*L)
			{
				tau0[k] = tauL;
			}
			if(x > x_c - 0.5*L && x <= x_c)
			{
				tau0[k] = tauL+(m-tauL)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			}
			if(x > x_c && x <= x_c + 0.5*L)
			{
				tau0[k] = tauR+(m-tauR)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			}
			if(x > x_c + 0.5*L)
			{
				tau0[k] = tauR;
			}
		}
		*/
		tauL = distribution_tauL(generator);
		tauR = distribution_tauR(generator);
		uR = distribution_uR(generator);
		uL = distribution_uL(generator);
		param[0][iter_MC] = tauL;
		param[1][iter_MC] = tauR;
		param[2][iter_MC] = uL;
		param[3][iter_MC] = uR;
		for (int k = 0; k < N; ++k)
		{
			double x = (k+0.5)*dx;
			if(x <= x_c)
			{
				tau0[k] = tauL;
				u0[k] = uL;
			}
			else
			{
				tau0[k] = tauR;
				u0[k] = uR;
			}
		}
		
		for (int i = 0; i < NP; ++i)
			param_out << param[i][iter_MC] << "\t";
		param_out << endl;
		
		roe_II st(tau0,u0,s_tau0,s_u0,gamma);
		st.set_CD(CD);
		TS.solve(st, true);
		st.get_U(W);
		for (int i = 0; i < D/2; ++i)
			for (int j = 0; j < N; ++j)
				average[i][j] += W[i][j];
		if (iter_MC < 1000)
		{
			st.print_physical(path+"samples_", ios::out | ios::app);
		}
	}
	
	for (int i = 0; i < D/2; ++i)
		for (int j = 0; j < N; ++j)
			average[i][j] /= MC_max;
	for (int j = 0; j < N; ++j)
	{
		tau_ave_out << average[0][j] << "\t";
		u_ave_out << average[1][j] << "\t";
	}

	for(int iter_MC = 0; iter_MC < MC_max; ++iter_MC)
	{
		cout << iter_MC << endl;
		/*
		tauL = param[0][iter_MC];
		tauR = param[1][iter_MC];
		m = param[2][iter_MC];
		L = param[3][iter_MC];
		x_c = param[4][iter_MC];
		
		for (int k = 0; k < N; ++k)
		{
			double x = (k+0.5)*dx;
			if(x <= x_c - 0.5*L)
			{
				tau0[0] = tauL;
			}
			if(x > x_c - 0.5*L && x <= x_c)
			{
				tau0[k] = tauL+(m-tauL)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			}
			if(x > x_c && x <= x_c + 0.5*L)
			{
				tau0[k] = tauR+(m-tauR)*sin(pi/L*(x-x_c)+pi/2)*sin(pi/L*(x-x_c)+pi/2);
			}
			if(x > x_c + 0.5*L)
			{
				tau0[k] = tauR;
			}
		}
		*/
		tauL = param[0][iter_MC];
		tauR = param[1][iter_MC];
		uL = param[2][iter_MC];
		uR = param[3][iter_MC];
		for (int k = 0; k < N; ++k)
		{
			double x = (k+0.5)*dx;
			if(x <= x_c)
			{
				tau0[k] = tauL;
				u0[k] = uL;
			}
			else
			{
				tau0[k] = tauR;
				u0[k] = uR;
			}
		}
		roe_II st(tau0,u0, s_tau0,s_u0, gamma);
		st.set_CD(CD);
		TS.solve(st, true);
		st.get_U(W);
		for (int i = 0; i < D/2; ++i)
			for (int j = 0; j < N; ++j)
				variance[i][j] += (W[i][j] - average[i][j])*(W[i][j] - average[i][j]);
	}
	
	for (int i = 0; i < D/2; ++i)
		for (int j = 0; j < N; ++j)
			variance[i][j] /= MC_max-1;
	for (int j = 0; j < N; ++j)
	{
		tau_var_out << variance[0][j] << "\t";
		u_var_out << variance[1][j] << "\t";
	}
	/***************************************/

	return 0;
}
