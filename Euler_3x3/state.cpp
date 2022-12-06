#include"state.hpp"

state::state(const vector<double>& rho, const vector<double>& u, const vector<double>& p, const vector<double>& s_rho, const vector<double>& s_u, const vector<double>& s_p, double g, const vector<double>& H, const vector<double>& delta_H, const vector<double>& S_H, const vector<double>& delta_S_H, bool c, bool d)
{
	D = 6;
	U.resize(D);
	for (int k = 0; k < D; ++k)
		U[k].resize(rho.size());

	U[0] = rho;
	U[3] = s_rho;
	gamma = g;
	h = H;
	delta_h = delta_H;
	s_h = S_H;
	delta_s_h = delta_S_H;
	for (unsigned int i = 0; i < rho.size(); ++i)
	{
		U[1][i] = rho[i]*u[i];
		U[2][i] = 0.5*rho[i]*u[i]*u[i] + p[i]/(gamma-1);
		U[4][i] = s_rho[i]*u[i] + rho[i]*s_u[i];
		U[5][i] = 0.5*s_rho[i]*u[i]*u[i] + rho[i]*u[i]*s_u[i] + s_p[i]/(gamma-1);
	}
	CD = c;
    activate_source_term = d;
	bc_L.assign(3,false);
	bc_R.assign(3,false);
};

state::state(const vector<double>& rho, const vector<double>& u, const vector<double>& p, const vector<double>& s_rho, const vector<double>& s_u, const vector<double>& s_p, double g, const vector<double>& H, const vector<double>& delta_H, bool c, bool d)
{
	D = 6;
	U.resize(D);
	for (int k = 0; k < D; ++k)
		U[k].resize(rho.size());
	
	U[0] = rho;
	U[3] = s_rho;
	gamma = g;
	h = H;
	delta_h = delta_H;
	s_h.assign(rho.size(), 0.);
	delta_s_h.assign(rho.size(), 0.);
	for (unsigned int i = 0; i < rho.size(); ++i)
	{
		U[1][i] = rho[i]*u[i];
		U[2][i] = 0.5*rho[i]*u[i]*u[i] + p[i]/(gamma-1);
		U[4][i] = s_rho[i]*u[i] + rho[i]*s_u[i];
		U[5][i] = 0.5*s_rho[i]*u[i]*u[i] + rho[i]*u[i]*s_u[i] + s_p[i]/(gamma-1);
	}
	CD = c;
    activate_source_term = d;
	bc_L.assign(3,false);
	bc_R.assign(3,false);
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
// get physical variables in cell k
void state::get_W(vector<double>& a, int k) const
{
	a.resize(D);
	a[0] = U[0][k];
	a[1] = U[1][k]/U[0][k]; // u = u*rho/rho
	a[2] = (gamma-1)*(U[2][k] - 0.5*U[1][k]*U[1][k]/U[0][k]); // p = (gamma-1)(E - 0.5(rho*u)^2/rho) = (gamma-1)(E - 0.5rho*u^2)
	a[3] = U[3][k];
	a[4] = (U[4][k] - U[3][k]*a[1])/U[0][k]; // s_u = (s_(rhou) - s_rho u)/rho
	a[5] = (gamma-1)*(U[5][k]-0.5*U[3][k]*a[1]*a[1] - U[1][k]*a[4]); //s_p = (gamma-1)(s_rhoE - 0.5*s_rho*u^2 - rho*u*s_u)
	
	return;
};

void state::physical2conservative(const vector<double>& W, vector<double>& U_output) const
{
	U_output.resize(6);
	U_output[0] = W[0];
	U_output[1] = W[0]*W[1];
	U_output[2] = 0.5*W[0]*W[1]*W[1] + W[2]/(gamma-1);
	U_output[3] = W[3];
	U_output[4] = W[3]*W[1] + W[0]*W[4];
	U_output[5] = 0.5*W[3]*W[1]*W[1] + W[0]*W[1]*W[4] + W[5]/(gamma-1);
	return;
};
void state::conservative2physical(const vector<double>& U_input, vector<double>& W) const
{
	W.resize(6);
	W[0] = U_input[0];
	W[1] = U_input[1]/U_input[0]; // u = u*rho/rho
	W[2] = (gamma-1)*(U_input[2] - 0.5*U_input[1]*U_input[1]/U_input[0]); // p = (gamma-1)(E - 0.5(rho*u)^2/rho) = (gamma-1)(E - 0.5rho*u^2)
	W[3] = U_input[3];
	W[4] = (U_input[4] - U_input[3]*W[1])/U_input[0];
	W[5] = (gamma-1)*(U_input[5]-0.5*U_input[3]*W[1]*W[1] - U_input[1]*W[4]);
	return;
};

double state::compute_H(const vector<double>& V) const
{
	return V[2]/V[0] + ((gamma-1)*(V[2] - 0.5*V[1]*V[1]/V[0]))/V[0];
};

double state::compute_s_H(const vector<double>& V) const
{
	double H = this->compute_H(V);
	double s_u = V[4]/V[0] - V[3]*V[1]/(V[0]*V[0]);
	double s_p = (gamma-1)*(V[5] - 0.5*V[3]*V[1]*V[1]/V[0]/V[0] - V[1]*s_u);
	return (V[5] + s_p)/V[0] - V[3]/V[0]*H;
};

// compute physical flux from conservative variables U
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
};

void state::print_conservative(const string& path, ios_base::openmode mode, int prec)
{
	ofstream file_rho (path+"rho.dat", mode);
	ofstream file_rhou (path+"rhou.dat", mode);
	ofstream file_rhoE (path+"rhoE.dat", mode);
	ofstream file_s_rho (path+"s_rho.dat", mode);
	ofstream file_s_rhou (path+"s_rhou.dat", mode);
	ofstream file_s_rhoE (path+"s_rhoE.dat", mode);
	
	file_rho.precision(prec);
	file_rhou.precision(prec);
	file_rhoE.precision(prec);
	file_s_rho.precision(prec);
	file_s_rhou.precision(prec);
	file_s_rhoE.precision(prec);
	
	for (unsigned int k = 0; k < U[0].size() - 1; ++k)
	{
		file_rho << U[0][k] << "\t";
		file_rhou << U[1][k] << "\t";
		file_rhoE << U[2][k] << "\t";
		file_s_rho << U[3][k] << "\t";
		file_s_rhou << U[4][k] << "\t";
		file_s_rhoE << U[5][k] << "\t";
	}
	file_rho << U[0][U[0].size() - 1] << endl;
	file_rhou << U[1][U[0].size() - 1] << endl;
	file_rhoE << U[2][U[0].size() - 1] << endl;
	file_s_rho << U[3][U[0].size() - 1] << endl;
	file_s_rhou << U[4][U[0].size() - 1] << endl;
	file_s_rhoE << U[5][U[0].size() - 1] << endl;
	
	return;
};
void state::print_physical(const string& path, ios_base::openmode mode, int prec)
{
	ofstream file_u (path+"u.dat", mode);
	ofstream file_rho (path+"rho.dat", mode);
	ofstream file_p (path+"p.dat", mode);
	ofstream file_s_u (path+"s_u.dat", mode);
	ofstream file_s_rho (path+"s_rho.dat", mode);
	ofstream file_s_p (path+"s_p.dat", mode);
	
	file_u.precision(prec);
	file_rho.precision(prec);
	file_p.precision(prec);
	file_s_u.precision(prec);
	file_s_rho.precision(prec);
	file_s_p.precision(prec);
	
	vector<vector<double> > W;
	this->get_W(W);
	
	for (unsigned int k = 0; k < W[0].size() - 1; ++k)
	{
		file_rho << W[0][k] << "\t";
		file_u << W[1][k] << "\t";
		file_p << W[2][k] << "\t";
		file_s_rho << W[3][k] << "\t";
		file_s_u << W[4][k] << "\t";
		file_s_p << W[5][k] << "\t";
	}
	file_rho << W[0][W[0].size() - 1] << endl;
	file_u << W[1][W[0].size() - 1] << endl;
	file_p << W[2][W[0].size() - 1] << endl;
	file_s_rho << W[3][W[0].size() - 1] << endl;
	file_s_u << W[4][W[0].size() - 1] << endl;
	file_s_p << W[5][W[0].size() - 1] << endl;

	return;
};
