#include"state.hpp"

state::state(vector<double> tau, vector<double> u, vector<double> s_tau, vector<double> s_u, double g)
{
	U.resize(4);
	U[0] = tau;
	U[1] = u;
	U[2] = s_tau;
	U[3] = s_u;
	gamma = g;
};
