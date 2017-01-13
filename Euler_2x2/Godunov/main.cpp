#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include"mesh.hpp"
#include"state.hpp"
#include"flux.hpp"

using namespace std;

int main()
{
	double xa(0), xb(1), dx(2e-3), T(0.1), t(0), dt;
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(1), uR(0), tauL(0.4), tauR(0.5), gamma(1.4), x_c(0.5); //uL(0), uR(2.246510973213954), tauL(0.2), tauR(1), gamma(1.4), x_c(0.5); //uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4), x_c(0.5);
	double sigma2 ((uR - uL)/(tauL - tauR));
	double cfl(0.5);
	vector<double> u0(N, uR);
	vector<double> tau0(N, tauR);
	// sensitivities
	vector<double> s_u0(N, 0);//(N, 1.183215956619923);
	vector<double> s_tau0(N,0);
	
	vector<int> dirac1(2), dirac2(2);
	vector<vector<double> > dirac;
	
	cout.precision(15);
	
	for (int k=0; k < N/2; ++k)
	{
		u0[k] = uL;
		tau0[k] = tauL;
		s_tau0[k] = 0;
		s_u0[k] = 1;//-9.351212140372281;
	}
	/*************** i.c. sol exact ***************/
	/*
	double u_star(1.415808559710058), tau_star(0.202240111571936), s_u_star(0.441585160424181), s_tau_star(-0.048004091679340);
	double sigma1 ((uL - u_star)/(tau_star - tauL));
	double sigma2 ((uR - u_star)/(tau_star - tauR));
	double t0 (0.005);
	for (int k = 0; k < N; ++k)
	{
	 	double x = (k+0.5)*dx;
	 	if (x - x_c < sigma1*t0)
	 	{
	  		u0[k] = uL;
	 		tau0[k] = tauL;
	 		s_u0[k] = 1;
	 		s_tau0[k] = 0;
	 	}
	 	if (x - x_c >= sigma1*t0 && x - x_c < sigma2*t0)
	 	{
	 		u0[k] = u_star;
	 		tau0[k] = tau_star;
	 		s_u0[k] = s_u_star;
	 		s_tau0[k] = s_tau_star;
	 	}
	 	if (x - x_c >= sigma2*t0)
	 	{
	 		u0[k] = uR;
	 		tau0[k] = tauR;
	 		s_u0[k] = 0;
	 		s_tau0[k] = 0;
	 	}
	}
	*/
	
	state st(tau0, u0, s_tau0, s_u0, gamma);
	vector<double> zeros(N+1, 0);
	vector<vector<double> > U, Uold, R, S;
	
	st.get_U(Uold);
	st.get_U(U);
	
	vector<double> sR, sL;
	flux fl;
	
	// printing on file initial data

	ofstream file_u ("../../../results/Euler_2x2_Godunov/shock_shock/u.dat");
	ofstream file_tau ("../../../results/Euler_2x2_Godunov/shock_shock/tau.dat");
	ofstream file_s_u ("../../../results/Euler_2x2_Godunov/shock_shock/s_u.dat");
	ofstream file_s_tau ("../../../results/Euler_2x2_Godunov/shock_shock/s_tau.dat");
	ofstream file_t ("../../../results/Euler_2x2_Godunov/shock_shock/t.dat");
	
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
		if(cont%100==0)
		{
			cout << "t = " << t << endl;
		}
		fl.residual(st, R, S);

		/*** dt computation ***/
		st.compute_lambdaR(sR);
		dt = dx*0.5*cfl/(*max_element(sR.begin(),sR.end()));
		if(first_time && (t+dt)>T)
		{
			dt = T-t;
			first_time = false;
		}
		/**********************/
		
		
		for (int i=0; i<N; ++i)
		{
			for (int k=0; k<4; ++k)
			{
				U[k][i] = Uold[k][i] + dt/dx*R[k][i];
			//	U[k+2][i] = Uold[k+2][i] + dt/dx*R[k+2][i] /*+ dt/dx*S[k][i]*/;
			}
		}
		
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