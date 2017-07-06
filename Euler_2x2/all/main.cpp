#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include<chrono>
#include"mesh.hpp"
#include"state.hpp"
#include"godunov.hpp"
#include"roe.hpp"
#include"roe_I.hpp"
#include"roe_II.hpp"
#include"utilities.hpp"

using namespace std;

int main()
{
	double xa(0), xb(1), dx(2e-3), T(0.03), t(0), dt;
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(0), uR(0), tauL(0.7), tauR(0.2), gamma(1.4), x_c(0.5);//uL(-1.563415104628313), uR(-3), tauL(0.2), tauR(0.5), gamma(1.4), x_c(0.5);
	double cfl(0.45);
	vector<double> u0(N, uR);
	vector<double> tau0(N, tauR);
	vector<double> s_u0(N, 0);
	vector<double> s_tau0(N,0);
	cout.precision(15);

	for (int k=0; k < N*x_c; ++k)
	{
		u0[k] = uL;
		tau0[k] = tauL;
		s_tau0[k] = 0;//1;
		s_u0[k] = 1;//-9.351212140372281;
	}
	
	vector<vector<double> > U, Uold, U_int, R;
	
	/***** CD *****/
	vector<double> x_bar(N+1,0), sigma(N+1,0);
	x_bar[0] = xa; x_bar[N] = xb;
	vector<vector<double> > U_bar(4, u0);
	/**************/
	
	/********** Restoring ***************/
/*	t = 0.028675308623962;
	ifstream if_tau ("tau_int.txt");
	ifstream if_u ("u_int.txt");
	ifstream if_s_tau ("s_tau_int.txt");
	ifstream if_s_u ("s_u_int.txt");
	double dummy; int k(0);
	while(if_tau >> dummy)
	{
		tau0[k] = dummy;
		++k;
	}
	k = 0;
	while(if_u >> dummy)
	{
		u0[k] = dummy;
		++k;
	}
	k = 0;
	while(if_s_tau >> dummy)
	{
		s_tau0[k] = dummy;
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
	
	roe_II st(tau0,u0,s_tau0,s_u0,gamma);
	bool time_secondorder (true);
	bool CD (false);
	st.set_CD(CD);
	vector<double> lambda;
	
	ofstream file_u ("results/u.dat");
	ofstream file_tau ("results/tau.dat");
	ofstream file_s_u ("results/s_u.dat");
	ofstream file_s_tau ("results/s_tau.dat");
	ofstream file_t ("results/t.dat");
	
	file_u.precision(15);
	file_tau.precision(15);
	file_s_u.precision(15);
	file_s_tau.precision(15);
	file_t.precision(15);
	
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
	
	st.get_U(Uold);
	st.get_U(U);
	st.get_U(U_int);
	
	bool first_time (true);
	int cont(0);
	// time loop
	
	auto timeSTART=chrono::system_clock::now();

	while (t < T)
	{
		++cont;
		if(cont%100==0)
		{
			cout << "t = " << t << endl;
		}
		
		/*** dt computation ***/
		st.compute_lambda(lambda);
		dt = dx*cfl/(*max_element(lambda.begin(),lambda.end()));
		if(first_time && (t+dt)>T)
		{
			dt = T-t;
			first_time = false;
		}
		/**********************/
		
		if (time_secondorder)
		{
			if (CD)
			{
				/** staggered grid definition **/
				sigma.assign(N+1,0);
				for (int i = 1; i < N; ++i)
				{
					if(U[1][i] - U[1][i-1] < -1e-3) // shock
					{
						if(U[0][i] - U[0][i-1] < -1e-3 ) //1-shock
							sigma[i] = -lambda[i];
						if(U[0][i] - U[0][i-1] > 1e-3) //2-shock
							sigma[i] = lambda[i];
					}
					x_bar[i] = dx*i + sigma[i]*0.5*dt;
				}
				/********************************/
				st.set_sigma(sigma);
				st.compute_residual(R);

				/********* compute U_int ********/
				for (int i=0; i<N; ++i)
				{
					double dxi = dx;
					double coeff = 0.0;
					if (fabs(sigma[i]) > 1e-8 || fabs(sigma[i+1]) > 1e-8)
					{
						dxi = (x_bar[i+1]-x_bar[i]);
						coeff = (dx-dxi)/dxi;
					}
					for (int k = 0; k < 4; ++k)
						U_int[k][i] = Uold[k][i] + coeff*Uold[k][i] + 0.5*dt/dxi*R[k][i];
				}
				/********************************/
				
				/********* compute U_bar ********/
				st.set_U(U_int);
				st.compute_residual(R);
				
				for (int i=0; i<N; ++i)
					x_bar[i] = dx*i + sigma[i]*dt;
				
				for (int i=0; i<N; ++i)
				{
					double dxi = dx;
					double coeff = 0.0;
					if (fabs(sigma[i]) > 1e-8 || fabs(sigma[i+1]) > 1e-8)
					{
						dxi = (x_bar[i+1]-x_bar[i]);
						coeff = (dx-dxi)/dxi;
					}
					for (int k = 0; k < 4; ++k)
						U_bar[k][i] = Uold[k][i] + coeff*U_int[k][i] + dt/dxi*R[k][i];
				}
				/********************************/

				/*********** sampling ***********/
				double an;
				can(cont, an);
				
				for (int i=0; i<N; ++i)
				{
					for (int k=0; k<4; ++k)
					{
						if (an < dt/dx*max(0.0, sigma[i]))
							U[k][i] = U_bar[k][i-1];

						if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i];
						
						if (an > 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i+1];

					}
				}
				/********************************/
			}
			else
			{
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
				{
					for (int k=0; k<4; ++k)
						U_int[k][i] = Uold[k][i] + 0.5*dt/dx*R[k][i];
					//cout << i << " " << R[0][i] << "\n";
				}
				//cout << "\nsecondo mezzo time step\n";
				st.set_U(U_int);
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
					for (int k=0; k<4; ++k)
						U[k][i] = Uold[k][i] + dt/dx*R[k][i];
			}
		}
		else
		{
			if(CD)
			{
				/** staggered grid definition **/
				sigma.assign(N+1,0);
				for (int i = 1; i < N; ++i)
				{
					if(U[1][i] - U[1][i-1] < -1e-3) // shock
					{
						if(U[0][i] - U[0][i-1] < -1e-3) //1-shock
							sigma[i] = -lambda[i];

						if(U[0][i] - U[0][i-1] > 1e-3) //2-shock
							sigma[i] = lambda[i];

					}
					x_bar[i] = dx*i + sigma[i]*dt;
				}
				/********************************/
				st.set_sigma(sigma);
				st.compute_residual(R);
				/********* compute U_bar ********/
				for (int i=0; i<N; ++i)
				{
					double dxi = (x_bar[i+1]-x_bar[i]);
					for (int k = 0; k < 4; ++k)
						U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
				}
				/********************************/
				
				/*********** sampling ***********/
				double an;
				can(cont, an);
				
				for (int i=0; i<N; ++i)
				{
					for (int k=0; k<4; ++k)
					{
						if (an < dt/dx*max(0.0, sigma[i]))
							U[k][i] = U_bar[k][i-1];

						if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i];

						if (an > 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i+1];

					}
				}
				/********************************/
			}
			else
			{
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
					for (int k=0; k<4; ++k)
						U[k][i] = Uold[k][i] + dt/dx*R[k][i];
			}
		}
		
		Uold = U;
		st.set_U(U);
		t += dt;
		
		if (cont%25 == 0 || !first_time)
		{
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
	}
	
	auto timeEND=chrono::system_clock::now();
	auto delta=chrono::duration_cast<chrono::milliseconds>(timeEND-timeSTART);
	int sec((delta.count())/1000);
	int msec((delta.count()-sec*1000));
	cout << sec << "s " << msec << endl;
	/* Non ho la soluzione esatta qui
	cout << endl << endl << "\\begin{tikzpicture}\n\\begin{axis}[title={$u(x,T)$}, legend pos=south west, legend style={fill=none}] \n\\addplot[red, thick] coordinates {";
	for (int i=0; i < N; ++i)
		cout << "(" << (0.5+i)*dx << ", " << << ") ";
	
	cout <<"}; \n\\addplot[blue] coordinates {"
	for (int i=0; i < N; ++i)
		cout << "(" << (0.5+i)*dx << ", " << << ") ";
	cout <<"}; \n\\legend{Exact, Numerical} \n\\end{axis} \n\\end{tikzpicture}\n";*/
	return 0;
}