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
	
	vector<double> U_R_cell(4), U_L_cell(4), U_star_cell(4), F_L_cell(4), F_R_cell(4); double lambda_cell, s_lambda_cell; int d1_cell(0), d2_cell(0);
	vector<double> UL(4), UR(4), U_L(4), U_R(4);
	vector<double> FL(4), FR(4), Ustar(4), Ustar_old(4), F_L(4), F_R(4);

	for (int i = 0; i <= N; ++i)
	{
		d1_cell = 0; d2_cell = 0;
		st.get_UL_extrapolated(UL, i);
		st.get_UR_extrapolated(UR, i);
		int left = max(i-1,0);
		int right = min(i, N-1);
		st.get_U(U_L, left);
		st.get_U(U_R, right);
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FL[2] = -UL[3];
		FL[3] = -gamma*UL[2]*pow(UL[0], -gamma-1);
		
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		FR[2] = -UR[3];
		FR[3] = -gamma*UR[2]*pow(UR[0], -gamma-1);
		
		F_L[0] = -U_L[1];
		F_L[1] = pow(U_L[0], -gamma);
		F_L[2] = -U_L[3];
		F_L[3] = -gamma*U_L[2]*pow(U_L[0], -gamma-1);
		
		F_R[0] = -U_R[1];
		F_R[1] = pow(U_R[0], -gamma);
		F_R[2] = -U_R[3];
		F_R[3] = -gamma*U_R[2]*pow(U_R[0], -gamma-1);

		for (int k = 0; k < 2; ++k)
		{
			Ustar[k] = 0.5*(UL[k] + UR[k]) - 0.5*(FR[k]-FL[k])/lambdaR[i];
			Ustar[k+2] = 0.5*(UL[k+2] + UR[k+2]) - 0.5*(FR[k+2]-FL[k+2])/lambdaR[i] + 0.5*s_lambdaR[i]/lambdaR[i]*(UR[k]-Ustar[k])*d2[i] + 0.5*s_lambdaR[i]/lambdaR[i]*(UL[k]-Ustar[k])*d1[i-1];
			Ustar_old[k] = 0.5*(U_L[k] + U_R[k]) - 0.5*(F_R[k]-F_L[k])/lambdaR[i];
		}
		
		for (int k=0; k<4; ++k)
		{
			if (i < N)
			{
				//R[k][i] += 0.5*(FL[k] + FR[k]) - 0.5*lambdaR[i]*(UR[k] - UL[k]);
				R[k][i] += lambdaR[i]*(Ustar[k]-UR[k]);
			}
			if (i > 0)
			{
				//R[k][i-1] -= 0.5*(FL[k] + FR[k]) - 0.5*lambdaR[i]*(UR[k] - UL[k]);
				R[k][i-1] += lambdaR[i]*(Ustar[k]-UL[k]);
			}
		}
		/*
		for (int k = 0; k < 2; ++k)
		{
	 	  	if(i < N)
			{
				S[k][i] += 0.5*s_lambdaR[i]*(U_R[k] - Ustar_old[k])*d2[i] + 0.5*(-s_lambdaR[i]*(Ustar_old[k] - U_L[k]))*d1[i-1];
				//S[k][i] += 0.5*s_lambdaR[i]*(UR[k] - Ustar[k])*d2[i] + 0.5*(-s_lambdaR[i]*(Ustar[k] - UL[k]))*d1[i-1];
			}
	 	  	if(i > 0)
			{
				S[k][i-1] += 0.5*(-s_lambdaR[i]*(Ustar_old[k] - U_L[k]))*d1[i-1] + 0.5*s_lambdaR[i]*(U_R[k] - Ustar_old[k])*d2[i];
				//S[k][i-1] += 0.5*(-s_lambdaR[i]*(Ustar[k] - UL[k]))*d1[i-1] + 0.5*s_lambdaR[i]*(UR[k] - Ustar[k])*d2[i];
			}
		}
		*/
		
		if(i > 0)
		{
			U_R_cell = UL;
			F_L_cell[0] = -U_L_cell[1];
			F_L_cell[1] = pow(U_L_cell[0], -gamma);
			F_L_cell[2] = -U_L_cell[3];
			F_L_cell[3] = -gamma*U_L_cell[2]*pow(U_L_cell[0], -gamma-1);
			F_R_cell[0] = -U_R_cell[1];
			F_R_cell[1] = pow(U_R_cell[0], -gamma);
			F_R_cell[2] = -U_R_cell[3];
			F_R_cell[3] = -gamma*U_R_cell[2]*pow(U_R_cell[0], -gamma-1);
			
			if(fabs(U_L_cell[0]-U_R_cell[0]) < 1e-15)
			{
				lambda_cell = sqrt( gamma*pow(U_L_cell[0], -gamma-1) );
				s_lambda_cell = -gamma*(gamma+1)*pow(U_L_cell[0], -gamma-2)*(U_R_cell[2]+U_L_cell[2])/(4*lambda_cell);
			}
			else
			{
				lambda_cell = sqrt( -(pow(U_L_cell[0], -gamma) - pow(U_R_cell[0], -gamma)) / (U_L_cell[0]-U_R_cell[0]) );
				s_lambda_cell = - (-gamma*U_L_cell[2]*pow(U_L_cell[0], -gamma-1) + gamma*U_R_cell[2]*pow(U_R_cell[0], -gamma-1))/(2*lambda_cell*(U_L_cell[0]-U_R_cell[0])) - lambda_cell*(U_L_cell[2]-U_R_cell[2])/(2*(U_L_cell[0]-U_R_cell[0]));
			}
			
			for (int k = 0; k < 2; ++k)
			{
				U_star_cell[k] = 0.5*(U_L_cell[k] + U_R_cell[k]) - 0.5*(F_R_cell[k]-F_L_cell[k])/lambda_cell;
			}
			
			if ((U_L_cell[0] - U_star_cell[0]) > 1e-3 && (U_L_cell[1] - U_star_cell[1]) > 1e-3 )
				d1_cell = 1;
			
			if ( (U_R_cell[0] - U_star_cell[0] > 1e-3) && (U_star_cell[1] - U_R_cell[1] > 1e-3))
				d2_cell = 1;
			
			for (int k = 0; k < 2; ++k)
			{
				U_star_cell[k+2] = 0.5*(U_L_cell[k+2] + U_R_cell[k+2]) - 0.5*(F_R_cell[k+2]-F_L_cell[k+2])/lambda_cell + 0.5*s_lambda_cell/lambda_cell*(U_R_cell[k]-U_star_cell[k])*d2_cell + 0.5*s_lambda_cell/lambda_cell*(U_L_cell[k]-U_star_cell[k])*d1_cell;
			}
			/*
			for (int k = 0; k < 2; ++k)
			{
				S[k][i-1] += 0.5*(-s_lambda_cell*(U_star_cell[k] - U_L_cell[k]))*d1_cell + 0.5*s_lambda_cell*(U_R_cell[k] - U_star_cell[k])*d2_cell;
			}
			*/
			for (int k = 0; k < 4; ++k)
				R[k][i-1] += lambda_cell*(2*U_star_cell[k] - U_L_cell[k] - U_R_cell[k]);

		}
		U_L_cell = UR;
		
	}
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
		st.get_UL_extrapolated(UL, i);
		st.get_UR_extrapolated(UR, i);
		
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
		st.get_UL_extrapolated(UL, i);
		st.get_UR_extrapolated(UR, i);
		
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