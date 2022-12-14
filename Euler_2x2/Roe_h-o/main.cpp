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
	double xa(0), xb(1), dx(1e-2), T(0.03), t(0), dt;
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(0), uR(0), tauL(0.7), tauR(0.2), gamma(1.4), x_c(0.5);//uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4), x_c(0.5); //uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4), x_c(0.5);
	double sigma2 ((uR - uL)/(tauL - tauR));
	double cfl(0.95);
	vector<double> u0(N, uR);
	vector<double> tau0(N, tauR);
	// sensitivities
	vector<double> s_u0(N, 0);
	vector<double> s_tau0(N,0);
	
	vector<int> dirac1(2), dirac2(2);
	vector<vector<double> > dirac;
	
	cout.precision(15);
	
	for (int k=0; k < N/2; ++k)
	{
		u0[k] = uL;
		tau0[k] = tauL;
		s_tau0[k] = 1;
		s_u0[k] = -9.351212140372281;
	}
	
	state st(tau0, u0, s_tau0, s_u0, gamma);
	state st_int(tau0, u0, s_tau0, s_u0, gamma);
	vector<double> zeros(N+1, 0);
	vector<vector<double> > U, Uold, U_int, R, S;
	
	st.get_U(Uold);
	st.get_U(U);
	st.get_U(U_int);
	
	vector<double> sR, sL;
	flux fl;
	
	// printing on file initial data
	
	ofstream file_u ("../../../results/Euler_2x2_Roe_h-o/2_shock/u.dat");
	ofstream file_tau ("../../../results/Euler_2x2_Roe_h-o/2_shock/tau.dat");
	ofstream file_s_u ("../../../results/Euler_2x2_Roe_h-o/2_shock/s_u.dat");
	ofstream file_s_tau ("../../../results/Euler_2x2_Roe_h-o/2_shock/s_tau.dat");
	ofstream file_t ("../../../results/Euler_2x2_Roe_h-o/2_shock/t.dat");
	
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
		dt = dx*cfl/(*max_element(sR.begin(),sR.end()));
		if(first_time && (t+dt)>T)
		{
			dt = T-t;
			first_time = false;
		}
		/**********************/
		
		
		for (int i=0; i<N; ++i)
		{
			for (int k=0; k<2; ++k)
			{
				U_int[k][i] = Uold[k][i] + 0.5*dt/dx*R[k][i];
				U_int[k+2][i] = Uold[k+2][i] + 0.5*dt/dx*R[k+2][i] + 0.5*dt/dx*S[k][i];
			}
		}
		
		st_int.set_U(U_int);
		fl.residual(st_int, R, S);
		
		for (int i=0; i<N; ++i)
		{
			for (int k=0; k<2; ++k)
			{
				U[k][i] = Uold[k][i] + dt/dx*R[k][i];
				U[k+2][i] = Uold[k+2][i] + dt/dx*R[k+2][i] + dt/dx*S[k][i];
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