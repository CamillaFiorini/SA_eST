#include"state.hpp"

state::state(vector<double> tau, vector<double> u, vector<double> s_tau, vector<double> s_u, double g)
{
	U.resize(4);
	int N = tau.size();
	for (int k = 0; k < 4; ++k)
	{
		U[k].resize(N);
	}
	U[0] = tau;
	U[1] = u;
	U[2] = s_tau;
	U[3] = s_u;
	gamma = g;
};

//get variables in cell i
void state::get_U(vector<double>& a, int i)
{
	a.resize(4);
	for (int k = 0; k < 4; ++k)
		a[k] = U[k][i];
	
	return;
};


void state::get_tau(vector<double>& a)
{
	a = U[0];
	return;
};
double state::get_tau(int i)
{
	return U[0][i];
};
void state::get_u(vector<double>& a)
{
	a = U[1];
	return;
};

double state::get_u(int i)
{
	return U[1][i];
};

void state::get_s_tau(vector<double>& a)
{
	a = U[2];
	return;
};
double state::get_s_tau(int i)
{
	return U[2][i];
};

void state::get_s_u(vector<double>& a)
{
	a = U[3];
	return;
};

double state::get_s_u(int i)
{
	return U[3][i];
};

// compute lambda_L
void state::compute_lambdaL(vector<double>& l)
{
	int N = U[0].size();
	l.resize(N+1);
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		double tauL = U[0][left];
		double tauR = U[0][right];
		if(fabs(tauL-tauR) < 1e-15)
			l[i] = - sqrt( gamma*pow(tauL, -gamma-1) );
		else
			l[i] = - sqrt( -(pow(tauL, -gamma) - pow(tauR, -gamma)) / (tauL-tauR) );
	}
	
	return;
};

void state::compute_lambdaR(vector<double>& l)
{
	int N = U[0].size();
	l.resize(N+1);
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		double tauL = U[0][left];
		double tauR = U[0][right];
		if(fabs(tauL-tauR) < 1e-15)
			l[i] = sqrt( gamma*pow(tauL, -gamma-1) );
		else
			l[i] = sqrt( -(pow(tauL, -gamma) - pow(tauR, -gamma)) / (tauL-tauR) );
	}
	return;
};

// compute d(lambda)/dp
void state::compute_s_lambdaL(vector<double>& sl)
{
	vector<double> l;
	this->compute_lambdaL(l);
	
	int N = U[0].size();
	sl.resize(N+1);
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		double tauL = U[0][left];
		double tauR = U[0][right];
		double s_tauL = U[2][left];
		double s_tauR = U[2][right];
		if(fabs(tauL-tauR) < 1e-15)
			sl[i] = -gamma*(gamma+1)*pow(tauL, -gamma-2)*(s_tauR+s_tauL)/(4*l[i]);
		else
			sl[i] = - (-gamma*s_tauL*pow(tauL, -gamma-1) + gamma*s_tauR*pow(tauR, -gamma-1))/(2*l[i]*(tauL-tauR)) - l[i]*(s_tauL-s_tauR)/(2*(tauL-tauR));
	}
	
	return;
};

void state::compute_s_lambdaR(vector<double>& sl)
{
	vector<double> l;
	this->compute_lambdaR(l);
	double dx = 2e-3;
	int N = U[0].size();
	sl.resize(N+1);
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		double tauL = U[0][left];
		double tauR = U[0][right];
		double s_tauL = U[2][left];
		double s_tauR = U[2][right];
		
		if(fabs(tauL-tauR) < 1e-15)
		{
			sl[i] = -gamma*(gamma+1)*pow(tauL, -gamma-2)*(s_tauR+s_tauL)/(4*l[i]);
		}
		else
		{
			sl[i] = - (-gamma*s_tauL*pow(tauL, -gamma-1) + gamma*s_tauR*pow(tauR, -gamma-1))/(2*l[i]*(tauL-tauR)) - l[i]*(s_tauL-s_tauR)/(2*(tauL-tauR));
		}
	}
	
	return;
};

void state::compute_Ustar(vector<vector<double> >& Ustar)
{
	int N = U[0].size();
	Ustar.resize(4);
	for (int k=0; k<4; ++k)
		Ustar[k].assign(N+1, 0);

	vector<double> lambda, s_lambda;
	vector<int> d1, d2;
	this->compute_lambdaR(lambda);
	this->compute_s_lambdaR(s_lambda);
	this->detector_s1(d1);
	this->detector_s2(d2);
	
	vector<double> FL(4), FR(4);
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		FL[0] = -U[1][left];
		FL[1] = pow(U[0][left], -gamma);
		FL[2] = -U[3][left];
		FL[3] = -gamma*U[2][left]*pow(U[0][left], -gamma-1);
	
		FR[0] = -U[1][right];
		FR[1] = pow(U[0][right], -gamma);
		FR[2] = -U[3][right];
		FR[3] = -gamma*U[2][right]*pow(U[0][right], -gamma-1);
		for (int k=0; k<2; ++k)
		{
			Ustar[k][i] = 0.5*(U[k][left] + U[k][right]) - 0.5/lambda[i]*(FR[k]-FL[k]);
			Ustar[k+2][i] = 0.5*(U[k+2][left] + U[k+2][right]) - 0.5/lambda[i]*(FR[k+2]-FL[k+2]) + 0.5*s_lambda[i]/lambda[i]*((U[k][right]-Ustar[k][i])*d2[right] + (U[k][left] - Ustar[k][i])*d1[left]);
		}
	}

	return;
};

void state::detector_s1(vector<int>& d1)
{
	int N = this->get_size();
	d1.assign(N, 0);
	vector<double> UL, UR, lambdaR, Ustar(2), FL(2), FR(2);
	this->compute_lambdaR(lambdaR);
	double gamma = this->get_gamma();
	
	for (int i = 0; i < N; ++i)
	{
		int left = i;
		int right = min(i+1, N-1);
		
		this->get_U(UL, left);
		this->get_U(UR, right);
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		
		for (int k=0; k<2; ++k)
			{
				Ustar[k] = 0.5*(UL[k]+UR[k]) - 0.5*(FR[k] - FL[k]) / lambdaR[i+1];
			}
		
		if ((UL[0] - Ustar[0]) > 0 && (UL[1] - Ustar[1]) > 0 )
			d1[i] = 1;
	}
	return;
}

void state::detector_s2(vector<int>& d2)
{
	int N = this->get_size();
	d2.assign(N, 0);
	vector<double> UL, UR, lambdaR, Ustar(2), FL(2), FR(2);
	this->compute_lambdaR(lambdaR);
	double gamma = this->get_gamma();
	
	for (int i = 0; i < N; ++i)
	{
		int left = max(i-1,0);
		int right = i;
		
		this->get_U(UL, left);
		this->get_U(UR, right);
		FL[0] = -UL[1];
		FL[1] = pow(UL[0], -gamma);
		FR[0] = -UR[1];
		FR[1] = pow(UR[0], -gamma);
		
		for (int k=0; k<2; ++k)
			{
				Ustar[k] = 0.5*(UL[k]+UR[k]) - 0.5*(FR[k] - FL[k]) / lambdaR[i];
			}
		
		
		if ( (UR[0] - Ustar[0] > 0) && (Ustar[1] - UR[1] > 0))
		{
			d2[i] = 1;
		}
	}
	return;
}