#include"state.hpp"

state::state(vector<double> tau, vector<double> u, vector<double> s_tau, vector<double> s_u, double g, bool c)
{
	U.resize(4);
	U[0] = tau;
	U[1] = u;
	U[2] = s_tau;
	U[3] = s_u;
	gamma = g;
    CD_state = c;
    CD_sens = c;
};

void state::set_U(vector<double> tau, vector<double> u, vector<double> s_tau, vector<double> s_u)
{
	U.resize(4);
	U[0] = tau;
	U[1] = u;
	U[2] = s_tau;
	U[3] = s_u;
};


void state::get_U(vector<double>& u, int i) const
{
	u.resize(4);
	for (int k = 0; k < 4; ++k)
		u[k] = U[k][i];
	return;
};

void state::print_physical(const string& path, ios_base::openmode mode, int prec)
{
	ofstream file_u (path+"u.dat", mode);
	ofstream file_tau (path+"tau.dat", mode);
	ofstream file_s_u (path+"s_u.dat", mode);
	ofstream file_s_tau (path+"s_tau.dat", mode);
	
	file_u.precision(prec);
	file_tau.precision(prec);
	file_s_u.precision(prec);
	file_s_tau.precision(prec);
	
	for (unsigned int k = 0; k < U[0].size(); ++k)
	{
		file_tau << U[0][k] << "\t";
		file_u << U[1][k] << "\t";
		file_s_tau << U[2][k] << "\t";
		file_s_u << U[3][k] << "\t";
	}
	file_tau << endl;
	file_u << endl;
	file_s_tau << endl;
	file_s_u << endl;
	
	return;
};
