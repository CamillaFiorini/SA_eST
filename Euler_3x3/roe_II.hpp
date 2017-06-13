#ifndef ROE_II_HPP
#define ROE_II_HPP

#include<vector>
#include"roe.hpp"
#include"utilities.hpp"

using namespace std;

class roe_II:public roe
{
protected:
	double kappa;
public:
	// Constructors
	roe_II(vector<vector<double> > u, double g) : roe(u,g), kappa(1./3.) {};
	roe_II(vector<double>a, vector<double>b, vector<double>c, vector<double>d, vector<double>e, vector<double>f, double g) : roe(a,b,c,d,e,f,g), kappa(1./3.) {};
	~roe_II() = default;
	// Methods
	void set_k(double k) {kappa=k;};
	double psi(double) const;
	void get_UL_extrapolated (vector<double>&, int) const;
	void get_UR_extrapolated (vector<double>&, int) const;
	void compute_residual(vector<vector<double> >&) const;
};
#endif