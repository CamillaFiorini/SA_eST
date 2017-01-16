#include"flux.hpp"

void flux::residual(state st, vector<vector<double> >& R, vector<vector<double> >& S)
{
	vector<vector<double> > F(4);
	vector<double> lambdaR, s_lambdaR;
	vector<int> d1, d2;
	st.compute_lambdaR(lambdaR);
	st.compute_s_lambdaR(s_lambdaR);
	double gamma = st.get_gamma();
	int N = st.get_size();
	R.resize(4);
	S.resize(2);
	S[0].assign(N,0);
	S[1].assign(N,0);
	for (int k=0; k<4; ++k)
	{
		F[k].resize(N+1);
		R[k].assign(N, 0);
	}
	
	this->detector_s1(st, d1, 1e-3);
	this->detector_s2(st, d2, 1e-3);
	
	cout.precision(15);
	
	for (int i = 0; i <= N; ++i)
	{
		int left = max(i-1,0);
		int Lleft = max(i-2, 0);
		int right = min(i, N-1);
		int Rright = min(i+1, N-1);
		
		vector<double> UL(4), UR(4), U_L, U_R, U_LL, U_RR;
		
		st.get_U(U_L, left);
		st.get_U(U_R, right);
		st.get_U(U_LL, Lleft);
		st.get_U(U_RR, Rright);
		double kappa = 1./3.;

		vector<double> rR(4), rL(4), irR(4), irL(4);
		
		for (int k = 0; k < 4; ++k)
		{
			if (fabs(U_RR[k]-U_R[k]) < 1e-15 || fabs((U_R[k]-U_L[k])) < 1e-15 /*|| d1[left] == 1 || d2[right] == 1*/)
			{
				irR[k] = 0;
				rR[k] = 0;
			}
			else
			{
				rR[k] = (U_R[k]-U_L[k])/(U_RR[k]-U_R[k]);
				irR[k] = 1./rR[k];
			}
			
			if (fabs(U_L[k]-U_LL[k]) < 1e-15 || fabs(U_R[k]-U_L[k]) < 1e-15 /*|| d1[left] == 1 || d2[right] == 1*/)
			{
				rL[k] = 0;
				irL[k] = 0;
			}
			else
			{
				rL[k] = (U_R[k]-U_L[k])/(U_L[k]-U_LL[k]);
				irL[k] = 1./rL[k];
			}
			UL[k] = U_L[k] +1.*( st.psi(rL[k])*(1-kappa)/4*(U_L[k] - U_LL[k]) + st.psi(irL[k])*(1+kappa)/4*(U_R[k] - U_L[k]));
			UR[k] = U_R[k] +1.*(- st.psi(irR[k])*(1+kappa)/4*(U_R[k] - U_L[k]) - st.psi(rR[k])*(1-kappa)/4*(U_RR[k] - U_R[k]));
			
		}
		
	//	if (i < N)
	//		cout << d2[i] << "\t";

	//	cout << st.psi(irR[0]) << "\t";
	/*	if(i > 20 && i < 30)
		{
			cout << "U_LL[0] = " << U_LL[0] << "\tU_L[0] = " << U_L[0] <<"\tU_R[0] = " << U_R[0] <<"\tU_RR[0] = " << U_RR[0] << "\trL[0] = " << rL[0] << "\trR[0] = " << rR[0] << "\t psi(rL[0]) = " << st.psi(rL[0]) << "\t psi(rR[0]) = " << st.psi(rR[0]) << endl;
		}
	*/	vector<double> FL(4), FR(4), Ustar(2); // 21.12 // Ustar(4);
		
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FL[2] = -UL[3];
		FL[3] = -gamma*UL[2]*pow(UL[0], -gamma-1);
		
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		FR[2] = -UR[3];
		FR[3] = -gamma*UR[2]*pow(UR[0], -gamma-1);
		
		for (int k = 0; k < 2; ++k)
		{
			Ustar[k] = 0.5*(UL[k] + UR[k]) - 0.5*(FR[k]-FL[k])/lambdaR[i];
			
			// 21.12 //
			//Ustar[k+2] = 0.5*(UL[k+2] + UR[k+2]) - 0.5*(FR[k+2]-FL[k+2])/lambdaR[i] + 0.5*s_lambdaR[i]/lambdaR[i]*((UR[k]-Ustar[k])*d2[right] + (UL[k]-Ustar[k])*d1[left]);
		}
		
		for (int k=0; k<4; ++k)
		{
			if (i < N)
			{
				R[k][i] += 0.5*(FL[k] + FR[k]) - 0.5*lambdaR[i]*(UR[k] - UL[k]);
				//R[k][i] += lambdaR[i]*(Ustar[k]-UR[k]);
			}
			if (i > 0)
			{
				R[k][i-1] -= 0.5*(FL[k] + FR[k]) - 0.5*lambdaR[i]*(UR[k] - UL[k]);
				//R[k][i-1] += lambdaR[i]*(Ustar[k]-UL[k]);
			}
		}
		
		for (int k = 0; k < 2; ++k)
		{
	 	  	if(i < N)
	 	  		S[k][i] += 0.5*s_lambdaR[i]*(UR[k] - Ustar[k])*d2[i] + 0.5*(-s_lambdaR[i]*(Ustar[k] - UL[k]))*d1[i-1];
	 	  	if(i > 0)
				S[k][i-1] += 0.5*(-s_lambdaR[i]*(Ustar[k] - UL[k]))*d1[i-1] + 0.5*s_lambdaR[i]*(UR[k] - Ustar[k])*d2[i];
		}
	}
//	cout << endl;
	return;
};

void flux::detector_s1(state st, vector<int>& d1, double threshold)
{
	int N (st.get_size());
	d1.assign(N, 0);
	vector<double> UL, UR, lambdaR, Ustar(2), FL(2), FR(2);
	st.compute_lambdaR(lambdaR);
	double gamma (st.get_gamma());
	
	for (int i = 0; i < N; ++i)
	{
		int left = max(i,0);
		int Lleft = max(i-1, 0);
		int right = min(i+1, N-1);
		int Rright = min(i+2, N-1);
		
		vector<double> UL(2), UR(2), U_L, U_R, U_LL, U_RR;
		
		st.get_U(U_L, left);
		st.get_U(U_R, right);
		st.get_U(U_LL, Lleft);
		st.get_U(U_RR, Rright);
		double kappa = 1./3.;
		
		vector<double> rR(2), rL(2), irR(2), irL(2);
		
		for (int k = 0; k < 2; ++k)
		{
			if (fabs(U_RR[k]-U_R[k]) < 1e-15 || fabs((U_R[k]-U_L[k])) < 1e-15)
			{
				irR[k] = 0;
				rR[k] = 0;
			}
			else
			{
				rR[k] = (U_R[k]-U_L[k])/(U_RR[k]-U_R[k]);
				irR[k] = 1./rR[k];
			}
			
			if (fabs(U_L[k]-U_LL[k]) < 1e-15 || fabs(U_R[k]-U_L[k]) < 1e-15)
			{
				rL[k] = 0;
				irL[k] = 0;
			}
			else
			{
				rL[k] = (U_R[k]-U_L[k])/(U_L[k]-U_LL[k]);
				irL[k] = 1./rL[k];
			}
			UL[k] = U_L[k] +1.*( st.psi(rL[k])*(1-kappa)/4*(U_L[k] - U_LL[k]) + st.psi(irL[k])*(1+kappa)/4*(U_R[k] - U_L[k]));
			UR[k] = U_R[k] +1.*(- st.psi(irR[k])*(1+kappa)/4*(U_R[k] - U_L[k]) - st.psi(rR[k])*(1-kappa)/4*(U_RR[k] - U_R[k]));
		}
		
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		
		for (int k=0; k<2; ++k)
		{
			Ustar[k] = 0.5*(UL[k]+UR[k]) - 0.5*(FR[k] - FL[k]) / lambdaR[i+1];
		}
		
		if ((UL[0] - Ustar[0]) > threshold && (UL[1] - Ustar[1]) > threshold )
			d1[i] = 1;
	}
	return;
}

void flux::detector_s2(state st, vector<int>& d2, double threshold)
{
	int N (st.get_size());
	d2.assign(N, 0);
	vector<double> UL, UR, lambdaR, Ustar(2), FL(2), FR(2);
	st.compute_lambdaR(lambdaR);
	double gamma (st.get_gamma());
	
	for (int i = 0; i < N; ++i)
	{
		int left = max(i-1,0);
		int Lleft = max(i-2, 0);
		int right = min(i, N-1);
		int Rright = min(i+1, N-1);
		
		vector<double> UL(2), UR(2), U_L, U_R, U_LL, U_RR;
		
		st.get_U(U_L, left);
		st.get_U(U_R, right);
		st.get_U(U_LL, Lleft);
		st.get_U(U_RR, Rright);
		double kappa = 1./3.;
		
		vector<double> rR(2), rL(2), irR(2), irL(2);
		
		for (int k = 0; k < 2; ++k)
		{
			if (fabs(U_RR[k]-U_R[k]) < 1e-15 || fabs((U_R[k]-U_L[k])) < 1e-15)
			{
				irR[k] = 0;
				rR[k] = 0;
			}
			else
			{
				rR[k] = (U_R[k]-U_L[k])/(U_RR[k]-U_R[k]);
				irR[k] = 1./rR[k];
			}
			
			if (fabs(U_L[k]-U_LL[k]) < 1e-15 || fabs(U_R[k]-U_L[k]) < 1e-15)
			{
				rL[k] = 0;
				irL[k] = 0;
			}
			else
			{
				rL[k] = (U_R[k]-U_L[k])/(U_L[k]-U_LL[k]);
				irL[k] = 1./rL[k];
			}
			UL[k] = U_L[k] +1.*( st.psi(rL[k])*(1-kappa)/4*(U_L[k] - U_LL[k]) + st.psi(irL[k])*(1+kappa)/4*(U_R[k] - U_L[k]));
			UR[k] = U_R[k] +1.*(- st.psi(irR[k])*(1+kappa)/4*(U_R[k] - U_L[k]) - st.psi(rR[k])*(1-kappa)/4*(U_RR[k] - U_R[k]));
		}
		
		st.get_U(UL, left);
		st.get_U(UR, right);
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		
		for (int k=0; k<2; ++k)
		{
			Ustar[k] = 0.5*(UL[k]+UR[k]) - 0.5*(FR[k] - FL[k]) / lambdaR[i];
		}

		if ( (UR[0] - Ustar[0] > threshold) && (Ustar[1] - UR[1] > threshold))
		{
			d2[i] = 1;
		}
	}
	return;
}