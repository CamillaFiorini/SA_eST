#include"flux.hpp"
/*
void flux::compute(state st, vector<vector<double> >& F)
{
	return;
};
*/
void flux::residual(state st, vector<vector<double> >& R, vector<vector<double> >& S)
{
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
		
		vector<double> FL(4), FR(4), Ustar(2); // 21.12 // Ustar(4);
		
		FL[0] = -uL;
		FL[1] = pow(tauL, -gamma);
		FL[2] = -s_uL;
		FL[3] = -gamma*s_tauL*pow(tauL, -gamma-1);
		
		FR[0] = -uR;
		FR[1] = pow(tauR, -gamma);
		FR[2] = -s_uR;
		FR[3] = -gamma*s_tauR*pow(tauR, -gamma-1);
		
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
			
			/*if(i < 253 && i > 247)
			{
				cout << "Ustar[" << k << "][" << i << "] = " << Ustar[k] << endl;
				cout << "Ustar[" << k+2 << "][" << i << "] = " << 0.5*(UL[k+2] + UR[k+2]) - 0.5*(FR[k+2]-FL[k+2])/lambdaR[i] + S[k][i]/lambdaR[i] << endl;
			}*/
			
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
}

void flux::source(state st, vector<vector<double> >& dirac)
{
	int N = st.get_size();
	vector<vector<double> > source, U; // source = dx*S !
	dirac.resize(2);
	source.resize(2);
	for (int k=0; k<2; ++k)
	{
		dirac[k].assign(N,0);
		source[k].assign(N+1,0);
	}
	vector<double> lambdaL, lambdaR, s_lambdaL, s_lambdaR;
	st.compute_lambdaL(lambdaL);
	st.compute_lambdaR(lambdaR);
	st.compute_s_lambdaL(s_lambdaL);
	st.compute_s_lambdaR(s_lambdaR);
	double gamma = st.get_gamma();
	st.get_U(U);
	cout.precision(15);

	
	vector<double> Ustar(4), UL, UR, FL(4), FR(4);
	
	for (int i = 0; i <= N; ++i)
	{
		int left = max(i-1,0);
		int right = min(i, N-1);
		
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
		
		FL[0] = -uL;
		FL[1] = pow(tauL, -gamma);
		FL[2] = -s_uL;
		FL[3] = -gamma*s_tauL*pow(tauL, -gamma-1);
		
		FR[0] = -uR;
		FR[1] = pow(tauR, -gamma);
		FR[2] = -s_uR;
		FR[3] = -gamma*s_tauR*pow(tauR, -gamma-1);
		
		if(fabs(lambdaR[i]) < 1e-10 && fabs(lambdaL[i]) < 1e-10)
			for (int k=0; k<4; ++k)
			{
				Ustar[k] = U[k][left];
			}
		else
			for (int k=0; k<4; ++k)
			{
				Ustar[k] = (lambdaR[i]*U[k][right] - lambdaL[i]*U[k][left])/(lambdaR[i]-lambdaL[i]) - (FR[k] - FL[k]) /(lambdaR[i]-lambdaL[i]);
			}
		
		for (int k = 0; k < 2; ++k)
		{
			source[k][i] = s_lambdaR[i]*(U[k][right] - Ustar[k])*(Ustar[0] < tauR)*(uR < Ustar[1]) - s_lambdaR[i]*(Ustar[k] - U[k][left])*(tauL > Ustar[0])*(Ustar[1] < uL); // dx*s[i-1/2]
		}
	}
	
	for (int i = 0; i < N; ++i)
	{
		for (int k = 0; k < 2; ++k)
		{
			dirac[k][i] = source[k][i] + source[k][i+1];
		}
	}
	
};