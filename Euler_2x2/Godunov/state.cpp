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

// compute lambda_R
void state::compute_lambdaR(vector<double>& l)
{
	int N = U[0].size();
	l.resize(N+1);
	vector<double> ustar;
	for (int i = 0; i <= N; ++i)
	{
		this->compute_Ustar(i, ustar);
		l[i] = sqrt( gamma*pow(ustar[0], -gamma-1) );
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
	vector<double> ustar;
	for (int i = 0; i <= N; ++i)
	{
		this->compute_Ustar(i, ustar);
		sl[i] = -0.5*sqrt(gamma)*(gamma+1)*pow(ustar[0], -0.5*(gamma+3));
	}
	
	return;
};


void state::compute_Ustar(int i, vector<double>& Ustar)
{
	int N = this->get_size();
	int left = max(0,i-1);
	int right = min(i, N-1);
	vector<double> UL, UR;
	
	this->get_U(UL, left);
	double tauL = UL[0];
	double uL = UL[1];
	double s_tauL = UL[2];
	double s_uL = UL[3];
	
	this->get_U(UR, right);
	double tauR = UR[0];
	double uR = UR[1];
	double s_tauR = UR[2];
	double s_uR = UR[3];
	
	vector<double> state_star(2);	
	this->newton (0.5*(tauL+tauR), 1e-15, 200, UL, UR, state_star);
	
	double s_tau_star = ( this->df2_dtauR(state_star[0], UR)*s_tauR + s_uR - this->df1_dtauL(state_star[0], UL)*s_tauL - s_uL )/(this->df1(state_star[0], UL) - this->df2(state_star[0], UR));
	double s_u_star = this->df1(state_star[0], UL)*s_tau_star + this->df1_dtauL(state_star[0], UL)*s_tauL + s_uL;
	
	if (isnan(s_tau_star) || isnan(s_u_star))
	{
		cout << "this->df1(state_star[0], UL) " << this->df1(state_star[0], UL) << endl;
	}
	
	Ustar.resize(4);
	Ustar[0] = state_star[0];
	Ustar[1] = state_star[1];
	Ustar[2] = s_tau_star;
	Ustar[3] = s_u_star;
	
	return;
};

void state::newton (double x, double toll, int iter_max, vector<double>& UL, vector<double>& UR, vector<double>& Ustar)
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

double state::f1(double x, vector<double> UL)
{
	double tauL = UL[0];
	double uL = UL[1];
	
	if (x < tauL - 1e-15)
		return uL - sqrt(-(pow(x,-gamma) - pow(tauL,-gamma))*(x-tauL));
	else
		return uL + 2*sqrt(gamma)/(1-gamma)*(pow(x, 0.5*(1-gamma)) - pow(tauL,0.5*(1-gamma)));
};

double state::f2(double x, vector<double> UR)
{
	double tauR = UR[0];
	double uR = UR[1];
	if (x < tauR - 1e-15)
		return uR + sqrt(-(pow(x,-gamma) - pow(tauR,-gamma))*(x-tauR));
	else
		return uR + 2*sqrt(gamma)/(1-gamma)*(pow(tauR, 0.5*(1-gamma)) - pow(x,0.5*(1-gamma)) );};

double state::df1(double x, vector<double> UL)
{
	double tauL = UL[0];
	double uL = UL[1];
	
	if(x < tauL - 1e-15)
		return (-gamma*pow(x,-gamma-1)*(x-tauL) + pow(x,-gamma)-pow(tauL,-gamma))/(2*sqrt(-(pow(x,-gamma) - pow(tauL,-gamma))*(x-tauL)));
	else
		return sqrt(gamma)*pow(x,-0.5*(gamma+1));
};

double state::df2(double x, vector<double> UR)
{
	double tauR = UR[0];
	double uR = UR[1];
	
	if(x < tauR - 1e-15)
		return -(-gamma*pow(x,-gamma-1)*(x-tauR) + pow(x,-gamma)-pow(tauR,-gamma))/(2*sqrt(-(pow(x,-gamma) - pow(tauR,-gamma))*(x-tauR)));
	else
		return -sqrt(gamma)*pow(x,-0.5*(gamma+1));
};

double state::df1_dtauL(double x, vector<double> UL)
{
	double tauL = UL[0];
	double uL = UL[1];
	
	if(x < tauL - 1e-15)
		return (gamma*pow(tauL,-gamma-1)*(x-tauL) - (pow(x,-gamma)-pow(tauL,-gamma)))/(2*sqrt(-(pow(x,-gamma) - pow(tauL,-gamma))*(x-tauL)));
	else
		return -sqrt(gamma)*pow(tauL,-0.5*(gamma+1));
};

double state::df2_dtauR(double x, vector<double> UR)
{
	double tauR = UR[0];
	double uR = UR[1];
	
	if(x < tauR - 1e-15)
		return (-gamma*pow(tauR,-gamma-1)*(x-tauR) + pow(x,-gamma)-pow(tauR,-gamma))/(2*sqrt(-(pow(x,-gamma) - pow(tauR,-gamma))*(x-tauR)));
	else
		return sqrt(gamma)*pow(tauR,-0.5*(gamma+1));
};