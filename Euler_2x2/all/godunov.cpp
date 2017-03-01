#include"godunov.hpp"
void godunov::compute_lambda(vector<double>& l) const
{
	int N = U[0].size();
	l.resize(N+1);
	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		vector<double> UL, UR, Ustar;
		this->get_U(UL, left);
		this->get_U(UR, right);
		this->compute_U_star(UL, UR, Ustar);
		l[i] = sqrt( gamma*pow(Ustar[0], -gamma-1) );
	}
	return;
};
void godunov::compute_U_star(const vector<double>& UL, const vector<double>& UR, vector<double>& Ustar) const
{
	double tauL = UL[0];
	double s_tauL = UL[2];
	double s_uL = UL[3];
	
	double tauR = UR[0];
	double s_tauR = UR[2];
	double s_uR = UR[3];
	
	vector<double> state_star(2);
	this->newton (0.5*(tauL+tauR), 1e-15, 200, UL, UR, state_star);
	
	double s_tau_star = ( this->df2_dtauR(state_star[0], UR)*s_tauR + s_uR - this->df1_dtauL(state_star[0], UL)*s_tauL - s_uL )/(this->df1(state_star[0], UL) - this->df2(state_star[0], UR));
	double s_u_star = this->df1(state_star[0], UL)*s_tau_star + this->df1_dtauL(state_star[0], UL)*s_tauL + s_uL;
	
	Ustar.resize(4);
	Ustar[0] = state_star[0];
	Ustar[1] = state_star[1];
	Ustar[2] = s_tau_star;
	Ustar[3] = s_u_star;
	
	return;
};

void godunov::compute_residual(vector<vector<double> >& R) const
{
	R.resize(4);
	int N = this->get_size();
	for (int k = 0; k < 4; ++k)
		R[k].assign(N,0);
	
 	for (int i = 0; i <= N; ++i)
	{
		int left = max(0,i-1);
		int right = min(i, N-1);
		vector<double> UL, UR, Ustar;
		this->get_U(UL, left);
		this->get_U(UR, right);
		this->compute_U_star(UL, UR, Ustar);
		
		double k1(0), k2(0);
		if(i > 0)
		{
			if( (Ustar[0]-UL[0] < -1e-15) && (Ustar[1] - UL[1] < -1e-15) ) // if 1-shock
				k1 = (Ustar[1] - UL[1])/(UL[0] - Ustar[0]);
			if ( (Ustar[0]-UL[0] > 1e-15) && (Ustar[1] - UL[1] > 1e-15) ) // if 1-raref
				k1 = -sqrt(gamma*pow(Ustar[0], -gamma-1));
		}
		if (i < N)
		{
			if( (Ustar[0]-UR[0] < -1e-15) && (Ustar[1] - UR[1] > 1e-15) ) // if 2-shock
				k2 = (Ustar[1] - UR[1])/(UR[0] - Ustar[0]);
			if( (Ustar[0]-UR[0] > 1e-15) && (Ustar[1] - UR[1] < -1e-15) ) // if 2-raref
				k2 = sqrt(gamma*pow(Ustar[0], -gamma-1));
		}
		if (i < N)
		{
			R[0][i] += -Ustar[1];
			R[1][i] += pow(Ustar[0], -gamma);
			R[2][i] += k2*(Ustar[2] - UR[2]);
			R[3][i] += k2*(Ustar[3] - UR[3]);
		}
		if (i > 0)
		{
			R[0][i-1] -= -Ustar[1];
			R[1][i-1] -= pow(Ustar[0], -gamma);
			R[2][i-1] -= k1*(Ustar[2] - UL[2]);
			R[3][i-1] -= k1*(Ustar[3] - UL[3]);
		}
	}
	
	return;
};

void godunov::newton (double x, double toll, int iter_max, const vector<double>& UL, const vector<double>& UR, vector<double>& Ustar) const
{
	double x_old = x + 1;
	int cont = 0;
	
	while (fabs(x-x_old) > toll && cont < iter_max)
	{
		++cont;
		x_old = x;
		x = x_old - (f1(x_old, UL) - f2(x_old, UR))/(df1(x_old, UL)-df2(x_old, UR));
	}
	if (fabs(x-x_old) > toll)
		cerr << "Newton did not converge\n";
	
	Ustar[0] = x;
	Ustar[1] = f1(x, UL);
	if (isnan(Ustar[0]) || isnan(Ustar[1]))
	{
		cerr << "There's a nan\n";
	}
	return;
};

double godunov::f1(double x, const vector<double>& UL) const
{
	double tauL = UL[0];
	double uL = UL[1];
	
	if (x < tauL - 1e-15)
		return uL - sqrt(-(pow(x,-gamma) - pow(tauL,-gamma))*(x-tauL));
	else
		return uL + 2*sqrt(gamma)/(1-gamma)*(pow(x, 0.5*(1-gamma)) - pow(tauL,0.5*(1-gamma)));
};

double godunov::f2(double x, const vector<double>& UR) const
{
	double tauR = UR[0];
	double uR = UR[1];
	if (x < tauR - 1e-15)
		return uR + sqrt(-(pow(x,-gamma) - pow(tauR,-gamma))*(x-tauR));
	else
		return uR + 2*sqrt(gamma)/(1-gamma)*(pow(tauR, 0.5*(1-gamma)) - pow(x,0.5*(1-gamma)) );};

double godunov::df1(double x, const vector<double>& UL) const
{
	double tauL = UL[0];
	
	if(x < tauL - 1e-15)
		return (-gamma*pow(x,-gamma-1)*(x-tauL) + pow(x,-gamma)-pow(tauL,-gamma))/(2*sqrt(-(pow(x,-gamma) - pow(tauL,-gamma))*(x-tauL)));
	else
		return sqrt(gamma)*pow(x,-0.5*(gamma+1));
};

double godunov::df2(double x, const vector<double>& UR) const
{
	double tauR = UR[0];
	
	if(x < tauR - 1e-15)
		return -(-gamma*pow(x,-gamma-1)*(x-tauR) + pow(x,-gamma)-pow(tauR,-gamma))/(2*sqrt(-(pow(x,-gamma) - pow(tauR,-gamma))*(x-tauR)));
	else
		return -sqrt(gamma)*pow(x,-0.5*(gamma+1));
};

double godunov::df1_dtauL(double x, const vector<double>& UL) const
{
	double tauL = UL[0];
	
	if(x < tauL - 1e-15)
		return (gamma*pow(tauL,-gamma-1)*(x-tauL) - (pow(x,-gamma)-pow(tauL,-gamma)))/(2*sqrt(-(pow(x,-gamma) - pow(tauL,-gamma))*(x-tauL)));
	else
		return -sqrt(gamma)*pow(tauL,-0.5*(gamma+1));
};

double godunov::df2_dtauR(double x, const vector<double>& UR) const
{
	double tauR = UR[0];
	
	if(x < tauR - 1e-15)
		return (-gamma*pow(tauR,-gamma-1)*(x-tauR) + pow(x,-gamma)-pow(tauR,-gamma))/(2*sqrt(-(pow(x,-gamma) - pow(tauR,-gamma))*(x-tauR)));
	else
		return sqrt(gamma)*pow(tauR,-0.5*(gamma+1));
};