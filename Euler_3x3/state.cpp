#include"state.hpp"

state::state(const vector<double>& rho, const vector<double>& u, const vector<double>& p, const vector<double>& s_rho, const vector<double>& s_u, const vector<double>& s_p, double g, bool c)
{
	D = 6;
	U.resize(D);
	for (int k = 0; k < D; ++k)
		U[k].resize(rho.size());

	U[0] = rho;
	U[3] = s_rho;
	gamma = g;
	for (unsigned int i = 0; i < rho.size(); ++i)
	{
		U[1][i] = rho[i]*u[i];
		U[2][i] = 0.5*rho[i]*u[i]*u[i] + p[i]/(gamma-1);
		U[4][i] = s_rho[i]*u[i] + rho[i]*s_u[i];
		U[5][i] = 0.5*s_rho[i]*u[i]*u[i] + rho[i]*u[i]*s_u[i] + s_p[i]/(gamma-1);
	}
	CD = c;
};

// get conservative variables
void state::get_U(vector<double>& u, int i) const
{
	u.resize(D);
	for (int k = 0; k < D; ++k)
		u[k] = U[k][i];
	return;
};

// get physical variables
void state::get_W(vector<vector<double> >& a) const
{
	a = U; //a[0] = U[0], a[3] = U[3] for the others:
	for (unsigned int k = 0; k < U[0].size(); ++k)
	{
		a[1][k] = U[1][k]/U[0][k]; // u = u*rho/rho
		a[2][k] = (gamma-1)*(U[2][k] - 0.5*U[1][k]*U[1][k]/U[0][k]); // p = (gamma-1)(E - 0.5(rho*u)^2/rho) = (gamma-1)(E - 0.5rho*u^2)
		a[4][k] = (U[4][k] - U[3][k]*a[1][k])/U[0][k]; // s_u =
		a[5][k] = (gamma-1)*(U[5][k]-0.5*U[3][k]*a[1][k]*a[1][k] - U[1][k]*a[4][k]); //s_p = (gamma-1)(s_E - 0.5*s_rho*u^2 - rho*u*s_u)
	}
	return;
};

void state::flux(const vector<double>& u, vector<double>& F) const
{
	F.resize(u.size());
	F[0] = u[1];
	F[1] = u[1]*u[1]/u[0] + (gamma-1)*(u[2] - 0.5*u[1]*u[1]/u[0]);
	F[2] = u[1]/u[0]*(u[2] + (gamma-1)*(u[2] - 0.5*u[1]*u[1]/u[0]));
	F[3] = u[4];
	double s_u = (u[4] - u[3]*u[1]/u[0])/u[0];
	F[4] = (u[3]*u[1]*u[1])/(u[0]*u[0]) + 2*u[1]*s_u + (gamma-1)*(u[5]-0.5*u[3]*u[1]/u[0]*u[1]/u[0] - u[1]*s_u); //s_rho*u*u + 2*rho*u*s_u + s_p
	F[5] = s_u*(u[2] + (gamma-1)*(u[2] - 0.5*u[1]*u[1]/u[0]) ) + u[1]/u[0]*(u[5] + (gamma-1)*(u[5]-0.5*u[3]*u[1]/u[0]*u[1]/u[0] - u[1]*s_u)); //s_uL*(EL+pL) + uL*(s_EL+s_pL);
}