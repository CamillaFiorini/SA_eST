#ifndef ROE_HPP
#define ROE_HPP

#include<vector>
#include"state.hpp"

using namespace std;

class roe:public state
{
public:
	// Constructors
	roe(vector<vector<double> > u, double g) : state(u,g) {};
	roe(vector<double>a, vector<double>b, vector<double>c, vector<double>d, double e) : state(a,b,c,d,e) {};
	// Methods
	void compute_lambda(vector<double>&) const;
	void compute_s_lambda(vector<double>&) const;
	void compute_U_star(const vector<double>&, const vector<double>&, vector<double>&) const;
	virtual void get_UL_extrapolated (vector<double>&, int) const = 0;
	virtual void get_UR_extrapolated (vector<double>&, int) const = 0;
	void compute_residual(vector<double>&) const;
};
#endif