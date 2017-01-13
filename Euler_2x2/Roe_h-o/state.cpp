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


// compute lambda_R
void state::compute_lambdaR(vector<double>& l)
{
	int N = U[0].size();
	l.resize(N+1);
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int Lleft = max(i-2,0);
		int right = min(i, N-1);
		int Rright = min(i+1, N-1);
		double tau_L = U[0][left];
		double tau_R = U[0][right];
		double tau_LL = U[0][Lleft];
		double tau_RR = U[0][Rright];
		double kappa = 1./3.;
		double irR, rR, irL, rL;
		//vector<int> d1, d2;
		
		if (fabs(tau_RR-tau_R) < 1e-15 || fabs((tau_R-tau_L)) < 1e-15 /*|| d1[left] == 1 || d2[right] == 1*/)
		{
			irR = 0;
			rR = 0;
		}
		else
		{
			rR = (tau_R-tau_L)/(tau_RR-tau_R);
			irR = 1./rR;
		}
		
		if (fabs(tau_L-tau_LL) < 1e-15 || fabs(tau_R-tau_L) < 1e-15 /*|| d1[left] == 1 || d2[right] == 1*/)
		{
			rL = 0;
			irL = 0;
		}
		else
		{
			rL = (tau_R-tau_L)/(tau_L-tau_LL);
			irL = 1./rL;
		}
		
		double tauL = tau_L +1.*( this->psi(rL)*(1-kappa)/4*(tau_L - tau_LL) + this->psi(irL)*(1+kappa)/4*(tau_R - tau_L));
		double tauR = tau_R +1.*(- this->psi(irR)*(1+kappa)/4*(tau_R - tau_L) - this->psi(rR)*(1-kappa)/4*(tau_RR - tau_R));
		
		if(fabs(tauL-tauR) < 1e-15)
			l[i] = sqrt( gamma*pow(tauL, -gamma-1) );
		else
			l[i] = sqrt( -(pow(tauL, -gamma) - pow(tauR, -gamma)) / (tauL-tauR) );
	}
	return;
};

// compute d(lambda)/dp
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