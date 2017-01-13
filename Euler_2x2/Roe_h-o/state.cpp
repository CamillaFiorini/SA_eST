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

double state::psi(double r)
{
	if(r < 1e-15)
		return 0;
	
	//return 1;
	//return 0.5*(1 + fabs(r)/r)*min(1., fabs(r));
	//return max(max(0.,min(2.,r)),min(1.,2*r));
	return (r*r+r)/(1+r*r);
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