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
double state::psi(double r)
{
	if(r < 1e-15)
		return 0;
	
	//return 1;
	return 0.5*(1 + fabs(r)/r)*min(1., fabs(r));
	//return max(max(0.,min(2.,r)),min(1.,2*r));
	//return (r*r+r)/(1+r*r);
};

// Functions that return the extrapolated values at the interface i
void state::get_UL_extrapolated (vector<double>& UL, int i, double kappa)
{
	UL.resize(4);
	int N = this->get_size();
	int left = max(i-1,0);
	int Lleft = max(i-2, 0);
	int right = min(i, N-1);
	int Rright = min(i+1, N-1);
	
	vector<double> rL(4), irL(4);
	
	for (int k = 0; k < 4; ++k)
	{
		if (fabs(U[k][left]-U[k][Lleft]) < 1e-15 || fabs(U[k][right]-U[k][left]) < 1e-15)
		{
			rL[k] = 0;
			irL[k] = 0;
		}
		else
		{
			rL[k] = (U[k][right]-U[k][left])/(U[k][left]-U[k][Lleft]);
			irL[k] = 1./rL[k];
		}
		UL[k] = U[k][left] +( this->psi(rL[k])*(1-kappa)/4*(U[k][left] - U[k][Lleft]) + this->psi(irL[k])*(1+kappa)/4*(U[k][right] - U[k][left]));
	}/*
	for (int k = 0; k < 4; ++k)
	{
		double a = U[k][right] - U[k][left];
		double b = U[k][left] - U[k][Lleft];
		if (a*b < 1e-10 )
			UL[k] = U[k][left];
		else
			UL[k] = U[k][left] + 0.5*fabs(a)/a*min(fabs(a), fabs(b));
	}*/
	return;
};

void state::get_UR_extrapolated (vector<double>& UR, int i, double kappa)
{
	UR.resize(4);
	int N = this->get_size();
	int left = max(i-1,0);
	int Lleft = max(i-2, 0);
	int right = min(i, N-1);
	int Rright = min(i+1, N-1);
	
	vector<double> rR(4), irR(4);
	
	for (int k = 0; k < 4; ++k)
	{
		if (fabs(U[k][Rright]-U[k][right]) < 1e-15 || fabs((U[k][right]-U[k][left])) < 1e-15)
		{
			irR[k] = 0;
			rR[k] = 0;
		}
		else
		{
			rR[k] = (U[k][right]-U[k][left])/(U[k][Rright]-U[k][right]);
			irR[k] = 1./rR[k];
		}
		UR[k] = U[k][right] +1.*(- this->psi(irR[k])*(1+kappa)/4*(U[k][right] - U[k][left]) - this->psi(rR[k])*(1-kappa)/4*(U[k][Rright] - U[k][right]));
	}
	/*
	for (int k = 0; k < 4; ++k)
	{
		double a = U[k][Rright] - U[k][right];
		double b = U[k][right] - U[k][left];
		if (a*b < 1e-10 )
			UR[k] = U[k][right];
		else
			UR[k] = U[k][right] - 0.5*fabs(a)/a*min(fabs(a), fabs(b));
	}*/
	return;
};

// compute lambda_R
void state::compute_lambdaR(vector<double>& l)
{
	int N = U[0].size();
	l.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		double tauL (UL[0]), tauR(UR[0]);
		
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
	int N = U[0].size();
	sl.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		
		double tauL (UL[0]), tauR(UR[0]), s_tauL (UL[2]), s_tauR(UR[2]);
		
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