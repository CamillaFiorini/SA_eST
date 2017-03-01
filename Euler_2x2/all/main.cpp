#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include"mesh.hpp"
#include"state.hpp"
#include"godunov.hpp"
#include"roe.hpp"
#include"roe_I.hpp"
#include"roe_II.hpp"

using namespace std;

int main()
{
	double xa(0), xb(1), dx(5e-3), T(0.08), t(0), dt;
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(0), uR(1), tauL(0.4), tauR(0.5), gamma(1.4);//uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4);
	double cfl(0.95);
	vector<double> u0(N, uR);
	vector<double> tau0(N, tauR);
	vector<double> s_u0(N, 0);
	vector<double> s_tau0(N,0);
	
	for (int k=0; k < N/2; ++k)
	{
		u0[k] = uL;
		tau0[k] = tauL;
		s_tau0[k] = 0;//1;
		s_u0[k] = 1;//-9.351212140372281;
	}
	
	vector<vector<double> > U, Uold, R, S;

	godunov st(tau0,u0,s_tau0,s_u0,gamma);
	roe_I st1(tau0,u0,s_tau0,s_u0,gamma);
	roe_II st2(tau0,u0,s_tau0,s_u0,gamma);
	
	st.compute_lambda(tau0);
	
	vector<double> sR, sL;
	
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
		
		/*** dt computation ***/
		st.compute_lambda(sR);
		dt = dx*cfl/(*max_element(sR.begin(),sR.end()));
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
			}
		}
		
		Uold = U;
		t += dt;
	}
	return 0;
}