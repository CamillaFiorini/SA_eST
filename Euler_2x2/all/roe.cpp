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
		l[i] = this->compute_lambda(UL, UR);
	}
	return;
};

double roe::compute_lambda(const vector<double>& UL, const vector<double>& UR) const
{
	if(fabs(UL[0]-UR[0]) < 1e-15)
		return sqrt( gamma*pow(UL[0], -gamma-1) );
	else
		return sqrt( -(pow(UL[0], -gamma) - pow(UR[0], -gamma)) / (UL[0]-UR[0]) );
};

double roe::compute_s_lambda(const vector<double>& UL, const vector<double>& UR) const
{
	double l = this->compute_lambda(UL,UR);
	
	if(fabs(UL[0]-UR[0]) < 1e-15)
		return -gamma*(gamma+1)*pow(UL[0], -gamma-2)*(UR[2]+UL[2])/(4*l);
	else
		return - (-gamma*UL[2]*pow(UL[0], -gamma-1) + gamma*UR[2]*pow(UR[0], -gamma-1))/(2*l*(UL[0]-UR[0])) - l*(UL[2]-UR[2])/(2*(UL[0]-UR[0]));
};

int roe::detector_s1(const vector<double>& UL, const vector<double>& UR, double threshold) const
{
	vector<double> Ustar;
	this->compute_state_star(UL, UR, Ustar);
	
	if ((UL[0] - Ustar[0]) > threshold && (UL[1] - Ustar[1]) > threshold )
		return 1;
	return 0;
};

int roe::detector_s2(const vector<double>& UL, const vector<double>& UR, double threshold) const
{
	vector<double> Ustar;
	this->compute_state_star(UL, UR, Ustar);
	
	if ( (UR[0] - Ustar[0] > threshold) && (Ustar[1] - UR[1] > threshold))
		return 1;
	return 0;
};

void roe::compute_state_star(const vector<double>& UL, const vector<double>& UR, vector<double>& Ustar) const
{
	vector<double> FL(2), FR(2);
	FL[0] = -UL[1];
	FL[1] = pow(UL[0], -gamma);
	
	FR[0] = -UR[1];
	FR[1] = pow(UR[0], -gamma);
	
	double lambda = this->compute_lambda(UL,UR);
	
	Ustar.resize(2);
	
	for (int k = 0; k < 2; ++k)
		Ustar[k] = 0.5*(UL[k] + UR[k]) - 0.5*(FR[k]-FL[k])/lambda;
	
	return;
};

void roe::compute_U_star(const vector<double>& UL, const vector<double>& UR, vector<double>& Ustar) const
{
	vector<double> FL(2), FR(2);
	FL[0] = -UL[3];
	FL[1] = -gamma*UL[2]*pow(UL[0], -gamma-1);
	
	FR[0] = -UR[3];
	FR[1] = -gamma*UR[2]*pow(UR[0], -gamma-1);
	
	double lambda = this->compute_lambda(UL,UR);
	double s_lambda = this->compute_s_lambda(UL, UR);
	
	int d1 = this->detector_s1(UL, UR);
	int d2 = this->detector_s2(UL, UR);
	
	this->compute_state_star(UL, UR, Ustar);
	
	Ustar.resize(4);
	for (int k = 0; k < 2; ++k)
		Ustar[k+2] = 0.5*(UL[k+2] + UR[k+2]) - 0.5*(FR[k]-FL[k])/lambda + 0.5*s_lambda/lambda*(UR[k]-Ustar[k])*d2 + 0.5*s_lambda/lambda*(UL[k]-Ustar[k])*d1;
	
	return;
};






