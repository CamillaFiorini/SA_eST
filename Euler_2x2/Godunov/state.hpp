#ifndef STATE_HPP
#define STATE_HPP

#include<vector>
#include<cmath>
#include<iostream>

using namespace std;

class state
{
private:
	vector<vector<double> > U;
	double gamma;
public:
	state();
	state(vector<double>,vector<double>,vector<double>,vector<double>, double);
	double compute_s(int i);
	inline int get_size() const {return U[0].size();};
	inline void get_U(vector<vector<double> >& a) {a = U;};
	inline void set_U(vector<vector<double> >& a) {U = a;};
	inline void set_gamma(double g) {gamma=g;};
	inline double get_gamma() {return gamma;};
	void get_U(vector<double>&, int);
	void get_tau(vector<double>&);
	double get_tau(int);
	void get_u(vector<double>&);
	double get_u(int);
	void get_s_tau(vector<double>&);
	double get_s_tau(int);
	void get_s_u(vector<double>&);
	double get_s_u(int);
	
	void compute_lambdaR(vector<double>&);
	void compute_s_lambdaR(vector<double>&);
	
	void compute_Ustar(int, vector<double>&);
	void newton(double, double, int, vector<double>&, vector<double>&, vector<double>&);
	double f1(double, vector<double>);
	double f2(double, vector<double>);
	double df1(double, vector<double>);
	double df2(double, vector<double>);
	double df1_dtauL(double, vector<double>);
	double df2_dtauR(double, vector<double>);
};

#endif