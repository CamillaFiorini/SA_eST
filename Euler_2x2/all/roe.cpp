#include"roe.hpp"

void roe::compute_lambda(vector<double>& l) const
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

void roe::compute_s_lambda(vector<double>& sl) const
{
	vector<double> l;
	this->compute_lambda(l);
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

void roe::compute_residual(vector<double>& R) const
{
	return;
};

void roe::compute_U_star(const vector<double>&, const vector<double>&, vector<double>&) const
{
	return;
};