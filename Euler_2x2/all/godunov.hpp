#ifndef GODUNOV_HPP
#define GODUNOV_HPP

#include<vector>
#include"state.hpp"

using namespace std;

class godunov:public state
{
public:
	// Constructors
	godunov(vector<vector<double> > u, double g) : state(u,g) {};
	godunov(vector<double>a, vector<double>b, vector<double>c, vector<double>d, double e) : state(a,b,c,d,e) {};
	// Methods
	void compute_lambda(vector<double>&) const;
	void compute_U_star(const vector<double>&, const vector<double>&, vector<double>&) const;
	void compute_residual(vector<vector<double> >&) const;
	
	//Functions that shouldn't be methods but I'm lazy
	void newton(double, double, int, const vector<double>&, const vector<double>&, vector<double>&) const;
	double f1(double, const vector<double>&) const;
	double f2(double, const vector<double>&) const;
	double df1(double, const vector<double>&) const;
	double df2(double, const vector<double>&) const;
	double df1_dtauL(double, const vector<double>&) const;
	double df2_dtauR(double, const vector<double>&) const;
};
#endif