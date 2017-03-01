#include"godunov.hpp"
void godunov::compute_lambda(vector<double>&) const
{
	return;
};
void godunov::compute_U_star(const vector<double>&, const vector<double>&, vector<double>&) const
{
	return;
};
void godunov::compute_residual(vector<double>&) const
{
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
	//cout << "f1(x_old, UL) = " << f1(x_old, UL) << "\tf2(x_old, UR) = " << f2(x_old, UR) << "\tUL[0] =  " << UR[1] << "\tUR[1] = " << UL[1] << endl;
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