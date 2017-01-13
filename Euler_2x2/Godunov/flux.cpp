#include"flux.hpp"
/*
void flux::compute(state st, vector<vector<double> >& F)
{
	return;
};
*/
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
	
	this->detector_s1(st, d1);
	this->detector_s2(st, d2);
	
	cout.precision(15);
	
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		vector<double> UL, UR;
		
		st.get_U(UL, left);
		double tauL = UL[0];
		double uL = UL[1];
		double s_tauL = UL[2];
		double s_uL = UL[3];
		
		st.get_U(UR, right);
		double tauR = UR[0];
		double uR = UR[1];
		double s_tauR = UR[2];
		double s_uR = UR[3];
		
		vector<double> Ustar(4);
		//st.newton (0.5*(tauL+tauR), 1e-15, 200, UL, UR, Ustar);
		st.compute_Ustar(i, Ustar);

		double s_tau_star = Ustar[2];// ( st.df2_dtauR(Ustar[0], UR)*s_tauR + s_uR - st.df1_dtauL(Ustar[0], UL)*s_tauL - s_uL )/(st.df1(Ustar[0], UL) - st.df2(Ustar[0], UR));
		double s_u_star = Ustar[3];//st.df1(Ustar[0], UL)*s_tau_star + st.df1_dtauL(Ustar[0], UL)*s_tauL + s_uL;
		
		
		if (i < 235 && i > 232)
		{
		//	cout << "s_tau_star[" << i << "] = " << s_tau_star << "\ts_u_star[" << i << "] = " << s_u_star << endl;
		}
		
		double sigma1(0), sigma2(0);
		if (i > 0)
		{
			if ( !(fabs(tauL - Ustar[0]) < 1e-15 &&  fabs(uL - Ustar[1]) < 1e-15) )
			{
				if (d1[i-1] == 1 && fabs(tauL - Ustar[0]) > 1e-15)
				{
					sigma1 = (Ustar[1] - uL)/(tauL - Ustar[0]);
				}
				else
				{
					sigma1 = -sqrt(gamma*pow(Ustar[0], -gamma-1));
				}
			}
		}
		if (i < N)
		{
			if ( !(fabs(tauR - Ustar[0]) < 1e-15 &&  fabs(uR - Ustar[1]) < 1e-15 ) ) // c'è un problema
			{
				if (d2[i] == 1 && fabs(tauR - Ustar[0]) > 1e-15)
				{
						sigma2 = (Ustar[1] - uR)/(tauR - Ustar[0]) ;
				}
				else
				{
					sigma2 = sqrt(gamma*pow(Ustar[0], -gamma-1));
				}
			}
		}
		
		if (i < 253 && i > 247)
		{
		//	cout << "sigma1 = " << sigma1 << " sigma2  = " << sigma2 << endl;
		}
		
		
		if (i < N)
		{
			R[0][i] += -Ustar[1];
			R[1][i] += pow(Ustar[0], -gamma);
			R[2][i] += -s_u_star;//sigma2*(s_tau_star - UR[2]);
			R[3][i] += -gamma*s_tau_star*pow(Ustar[0], -gamma-1); //sigma2*(s_u_star - UR[3]);
		}
		if (i > 0)
		{
			R[0][i-1] -= -Ustar[1];
			R[1][i-1] -= pow(Ustar[0], -gamma);
			R[2][i-1] -= -s_u_star; //sigma1*(s_tau_star - UL[2]);
			R[3][i-1] -= -gamma*s_tau_star*pow(Ustar[0], -gamma-1); //sigma1*(s_u_star - UL[3]);
		}
		
		if (isnan(R[3][i]))
		{
			//cout << "s_u_star = " << s_u_star << endl;
		}
		
		for (int k = 0; k < 2; ++k)
		{
	 	  	if(i < N)
	 	  		S[k][i] += s_lambdaR[i]*(UR[k] - Ustar[k])*d2[i];
	 	  	if(i > 0)
	 	  		S[k][i-1] += -s_lambdaR[i]*(Ustar[k] - UL[k])*d1[i-1];
		}
	}
	return;
};

void flux::detector_s1(state st, vector<int>& d1)
{
	int N (st.get_size());
	d1.assign(N, 0);
	vector<double> UL, UR, lambdaR, Ustar(2), FL(2), FR(2);
	st.compute_lambdaR(lambdaR);
	double gamma (st.get_gamma());
	
	/********** Analytical **************/
//	double sigma1 = -3.239309819627467;
//	double dx = 2e-3;
	/************************************/

	
	for (int i = 0; i < N; ++i)
	{
		int left = i;
		int right = min(i+1, N-1);
		
		st.get_U(UL, left);
		st.get_U(UR, right);
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		
/*		if(fabs(lambdaR[i+1]) < 1e-15)
			for (int k=0; k<2; ++k)
			{
				Ustar[k] = UL[k];
			}
		else
*/			for (int k=0; k<2; ++k)
			{
				Ustar[k] = 0.5*(UL[k]+UR[k]) - 0.5*(FR[k] - FL[k]) / lambdaR[i+1];
			}
		
		if ((UL[0] - Ustar[0]) > 0 && (UL[1] - Ustar[1]) > 0 )
			d1[i] = 1;
 
		/********** Analytical **************/
/*		double xs = 0.5 + sigma1*t;
		if (xs >= dx*i && xs < dx*(i+1))
		{
			d1[i] = 1;
			d1[i-1] = 1;
			d1[i+1] = 1;
			d1[i-2] = 1;
			d1[i+2] = 1;
			d1[i-3] = 1;
			d1[i+3] = 1;
		}
*/		/************************************/
		
	}
	return;
}

void flux::detector_s2(state st, vector<int>& d2)
{
	
	/********** Analytical **************/
//	double sigma2 = 4.754866638298402;
//	double dx = 2e-3;
	/************************************/

	
	int N (st.get_size());
	d2.assign(N, 0);
	vector<double> UL, UR, lambdaR, Ustar(2), FL(2), FR(2);
	st.compute_lambdaR(lambdaR);
	double gamma (st.get_gamma());
	
	for (int i = 0; i < N; ++i)
	{
		int left = max(i-1,0);
		int right = i;
		
		st.get_U(UL, left);
		st.get_U(UR, right);
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		
/*		if(fabs(lambdaR[i]) < 1e-15)
			for (int k=0; k<2; ++k)
			{
				Ustar[k] = UL[k];
			}
		else
*/			for (int k=0; k<2; ++k)
			{
				Ustar[k] = 0.5*(UL[k]+UR[k]) - 0.5*(FR[k] - FL[k]) / lambdaR[i];
			}
		

		if ( (UR[0] - Ustar[0] > 0) && (Ustar[1] - UR[1] > 0))
		{
			d2[i] = 1;
		}
		
		/********** Analytical **************/
/*		double xs = 0.5 + sigma2*t;
		if (xs >= dx*i && xs < dx*(i+1))
		{
			d2[i] = 1;
			d2[i-1] = 1;
			d2[i+1] = 1;
			d2[i-2] = 1;
			d2[i+2] = 1;
			d2[i-3] = 1;
			d2[i+3] = 1;
		}
*/		/************************************/
		
		
	}
	return;
};
