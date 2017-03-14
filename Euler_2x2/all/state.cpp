#include"state.hpp"

state::state(vector<double> tau, vector<double> u, vector<double> s_tau, vector<double> s_u, double g, bool c)
{
	U.resize(4);
	U[0] = tau;
	U[1] = u;
	U[2] = s_tau;
	U[3] = s_u;
	gamma = g;
	CD = c;
};

void state::get_U(vector<double>& u, int i) const
{
	u.resize(4);
	for (int k = 0; k < 4; ++k)
		u[k] = U[k][i];
	return;
};
