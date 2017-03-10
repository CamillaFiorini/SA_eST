#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include"mesh.hpp"
#include"state.hpp"
#include"flux.hpp"
#include"utilities.hpp"

using namespace std;

int main()
{
	double xa(0), xb(1), dx(1e-2), T(0.1), t(0), dt;
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(0), uR(2.246510973213954), tauL(0.2), tauR(1), gamma(1.4), x_c(0.5);//uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4), x_c(0.5);//uL(0), uR(0), tauL(0.7), tauR(0.2), gamma(1.4), x_c(0.5);
	double cfl(0.95);
	// state vectors
	vector<double> u0(N, uR);
	vector<double> tau0(N, tauR);
	// sensitivities
	vector<double> s_u0(N, 1.183215956619923);//(N, 0);
	vector<double> s_tau0(N,1);
	
	/*****************************************************************************/
	/********************************** CD ***************************************/
	/*****************************************************************************/
	vector<double> x_bar(N+1,0), sigma(N+1,0);
	x_bar[0] = xa; x_bar[N] = xb;
	vector<vector<double> > U_bar(4, u0), F_bar(4, sigma);
	/*****************************************************************************/
	
	for (int k=0; k < N/2; ++k)
	{
		u0[k] = uL;
		tau0[k] = tauL;
		s_tau0[k] = 0;//1;
		s_u0[k] = 0;//-9.351212140372281;
	}
	
	state st(tau0, u0, s_tau0, s_u0, gamma);
	vector<double> zeros(N+1, 0);
	vector<vector<double> > U, Uold, R, S;
	
	st.get_U(Uold);
	st.get_U(U);
	
	vector<double> lambda;
	flux fl;
	
	// printing on file initial data
	ofstream file_u ("../../../results/Euler_2x2_CD_Roe_h-o/1_raref/u.dat");
	ofstream file_tau ("../../../results/Euler_2x2_CD_Roe_h-o/1_raref/tau.dat");
	ofstream file_s_u ("../../../results/Euler_2x2_CD_Roe_h-o/1_raref/s_u.dat");
	ofstream file_s_tau ("../../../results/Euler_2x2_CD_Roe_h-o/1_raref/s_tau.dat");
	ofstream file_t ("../../../results/Euler_2x2_CD_Roe_h-o/1_raref/t.dat");
	
	file_u.precision(15);
	file_tau.precision(15);
	file_s_u.precision(15);
	file_s_tau.precision(15);

	for (int k = 0; k < N; ++k)
	{
		file_u << u0[k] << "\t";
		file_tau << tau0[k] << "\t";
		file_s_u << s_u0[k] << "\t";
		file_s_tau << s_tau0[k] << "\t";
	}
	file_u << endl;
	file_tau << endl;
	file_s_tau << endl;
	file_s_u << endl;
	file_t << t << endl;
	
	bool first_time (true);
	int cont(0);
	// time loop
	while (t < T)
	{
		++cont;
		if(cont%20==0)
		{
			cout << "t = " << t << endl;
		}
		
		/********************************************************************************************************/
		/******************************************** dt computation ********************************************/
		/********************************************************************************************************/
		st.compute_lambdaR(lambda);
		dt = dx*cfl/(*max_element(lambda.begin(),lambda.end()));
		if(first_time && (t+dt)>T)
		{
			dt = T-t;
			first_time = false;
		}
		/********************************************************************************************************/
		
		/********************************************************************************************************/
		/****************************************** staggered grid def ******************************************/
		/********************************************************************************************************/
		vector<vector<double> > Ustar;
		vector<double> s_lambda;
		st.compute_s_lambdaR(s_lambda);
		st.compute_Ustar(Ustar);
		fl.residual(st, R, S);
		
		sigma.assign(N+1,0);
		for (int i = 1; i < N; ++i)
		{
			if(U[1][i] < U[1][i-1]) // shock
			{
				if(U[0][i] < U[0][i-1]) //1-shock
				{
					sigma[i] = -lambda[i];
				}
				if(U[0][i] > U[0][i-1]) //2-shock
				{
					sigma[i] = lambda[i];
				}
			}
			x_bar[i] = dx*i + sigma[i]*dt;
		}
		/********************************************************************************************************/

		/********************************************************************************************************/
		/******************************************** compute U_bar *********************************************/
		/********************************************************************************************************/
		for (int i=0; i<N; ++i)
		{
			double dxi = (x_bar[i+1]-x_bar[i]);
			for (int k = 0; k < 4; ++k)
			{
				U_bar[k][i] = 1.0/dxi*( (dx + dt*(-lambda[i+1]-lambda[i]))*Uold[k][i] + (sigma[i+1]+lambda[i+1])*dt*Ustar[k][i+1] + (lambda[i]-sigma[i])*dt*Ustar[k][i]);
				U[k][i] = Uold[k][i] + dt/dx*R[k][i];
			}
		}
		/********************************************************************************************************/

		/***********************/
		/****** compute an *****/
		/***********************/
		double an;
		can(cont, an);
		/***********************/
		/*
		for (int i = 0; i < N; ++i)
		{
			for (int k = 0; k < 4; ++k)
			{
				if (an < dt/dx*max(0.0, sigma[i]))
				{
					U[k][i] = U_bar[k][i-1];
				}
				if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
				{
					U[k][i] = U_bar[k][i];
				}
				if (an > 1+dt/dx*min(0.0, sigma[i+1]))
				{
					U[k][i] = U_bar[k][i+1];
				}
			}
		}
		*/
		Uold = U;
		st.set_U(U);
		t += dt;

		for (int k=0; k<N; ++k)
		{
			file_tau << U[0][k] << "\t";
			file_u << U[1][k] << "\t";
			file_s_tau << U[2][k] << "\t";
			file_s_u << U[3][k] << "\t";

		}
		file_tau << endl;
		file_u << endl;
		file_s_tau << endl;
		file_s_u << endl;
		file_t << t << endl;
	}
	return 0;
}