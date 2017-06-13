#ifndef ROE_HPP
#define ROE_HPP

#include<vector>
#include<algorithm>
#include"state.hpp"

using namespace std;

class roe:public state
{
public:
	// Constructors
	roe(vector<vector<double> > u, double g) : state(u,g) {};
	roe(vector<double>a, vector<double>b, vector<double>c, vector<double>d, vector<double>e, vector<double>f, double g) : state(a,b,c,d,e,f,g) {};
	// Methods
	double compute_lambda1(const vector<double>&, const vector<double>&) const;
	void compute_lambda1(vector<double>&) const;
	double compute_lambda2(const vector<double>&, const vector<double>&) const;
	void compute_lambda2(vector<double>&) const;
	double compute_lambda3(const vector<double>&, const vector<double>&) const;
	void compute_lambda3(vector<double>&) const;
	double compute_maxvel() const;
	double compute_utilde(const vector<double>&, const vector<double>&) const;
	double compute_s_utilde(const vector<double>&, const vector<double>&) const;
	double compute_H(const vector<double>&) const;
	double compute_s_H(const vector<double>&) const;
	double compute_Htilde(const vector<double>& UL, const vector<double>& UR) const;
	double compute_s_Htilde(const vector<double>& UL, const vector<double>& UR) const;
	double compute_atilde(const vector<double>&, const vector<double>&) const;
	double compute_s_atilde(const vector<double>&, const vector<double>&) const;
	void compute_alpha_tilde(const vector<double>&, const vector<double>&, vector<double>&) const;
	int detector_s1(const vector<double>&, const vector<double>&, double=1e-10) const;
	int detector_s2(const vector<double>&, const vector<double>&, double=1e-10) const;
	int detector_c(const vector<double>&, const vector<double>&, double=1e-10) const;
	void compute_flux(const vector<double>&, const vector<double>&, vector<double>&, vector<double>&) const;
	virtual void get_UL_extrapolated (vector<double>&, int) const = 0;
	virtual void get_UR_extrapolated (vector<double>&, int) const = 0;
	void compute_residual(vector<vector<double> >&) const;
};
#endif