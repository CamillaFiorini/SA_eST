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

using namespace std;

int main()
{
	double xa(0), xb(1), dx(1e-3), T(0.1), t(0), dt;
	mesh M (xa, xb, dx);
	int N = M.get_N();
	double uL(0), uR(0), rhoL(1), rhoR(0.125), pR(0.1), pL(1), x_c(0.5), gamma(1.4);//uL(19.5975), uR(-6.19633), rhoL(5.99924), rhoR(5.99242), pR(46.0950), pL(460.894), x_c(0.5), gamma(1.4);
	double cfl(0.5);
	vector<double> u0(N, uR);
	vector<double> rho0(N, rhoR);
	vector<double> p0(N, pR);
	vector<double> s_u0(N, 0);
	vector<double> s_rho0(N,1);
	vector<double> s_p0(N,0);
	cout.precision(15);

	for (int k=0; k < N*x_c; ++k)
	{
		u0[k] = uL;
		rho0[k] = rhoL;
		p0[k] = pL;
		s_rho0[k] = 0;
		s_u0[k] = 0;
		s_p0[k] = 0;
	}
	
	vector<vector<double> > U, Uold, U_int, W, R;
	
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
	roe_I st(rho0,u0, p0, s_rho0,s_u0, s_p0,gamma);
	bool time_secondorder (false);
	bool CD (false);
	st.set_CD(CD);
	
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
	
	st.get_U(Uold);
	st.get_U(U);
	st.get_U(U_int);
	
	bool first_time (true);
	int cont(0);
	// time loop
	
	auto timeSTART=chrono::system_clock::now();

	while (t < T)
	{
		vector<double> UL, UR;
		for (int i = 0; i < N+1; ++i)
		{
			st.get_UR_extrapolated(UR, i);
			st.get_UL_extrapolated(UL, i);
			file_d1 << st.detector_s1(UL, UR, 1e-3) << "\t";
			file_d2 << st.detector_s2(UL, UR, 1e-3) << "\t";
			file_c << st.detector_c(UL, UR, 1e-3) << "\t";
		}
		file_d1 << endl;
		file_d2 << endl;
		file_c << endl;
		++cont;
		if(cont%100==0)
		{
			cout << "t = " << t << endl;
		}
		
		/*** dt computation ***/
		double mv = st.compute_maxvel();
		dt = dx*cfl/mv;
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
	/*			sigma.assign(N+1,0);
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
	*/			/********************************/
	/*			st.set_sigma(sigma);
				st.compute_residual(R);
	*/
				/********* compute U_int ********/
	/*			for (int i=0; i<N; ++i)
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
	*/			/********************************/
				
				/********* compute U_bar ********/
	/*			st.set_U(U_int);
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
	*/			/********************************/

				/*********** sampling ***********/
	/*			double an;
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
	*/			/********************************/
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
	/*			sigma.assign(N+1,0);
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
	*/			/********************************/
	/*			st.set_sigma(sigma);
				st.compute_residual(R);
	*/			/********* compute U_bar ********/
	/*			for (int i=0; i<N; ++i)
				{
					double dxi = (x_bar[i+1]-x_bar[i]);
					for (int k = 0; k < 4; ++k)
						U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
				}
	*/			/********************************/
				
				/*********** sampling ***********/
	/*			double an;
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
	*/			}
				/********************************/
	//		}
			else
			{
				st.compute_residual(R);				
				for (int i=0; i<N; ++i)
					for (int k=0; k<6; ++k)
						U[k][i] = Uold[k][i] + dt/dx*R[k][i];
			}
		}
		
		Uold = U;
		st.set_U(U);
		t += dt;
		
		st.get_W(W);
		if (cont%10 == 0 || !first_time)
		{
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