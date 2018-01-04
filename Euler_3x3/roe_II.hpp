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
	roe_II(const vector<vector<double> >& u, double g) : roe(u,g), kappa(1./3.) {};
	roe_II(const vector<double>&a, const vector<double>& b, const vector<double>& c, const vector<double>& d, const vector<double>& e, const vector<double>& f, double g,const vector<double>& h, const vector<double>& dh) : roe(a,b,c,d,e,f,g,h,dh), kappa(1./3.) {};
	roe_II(const vector<double>&a, const vector<double>& b, const vector<double>& c, const vector<double>& d, const vector<double>& e, const vector<double>& f, double g,const vector<double>& h, const vector<double>& dh, const vector<double>& sh, const vector<double>& dsh) : roe(a,b,c,d,e,f,g,h,dh,sh,dsh), kappa(1./3.) {};
	~roe_II() = default;
	// Methods
	void set_k(double k) {kappa=k;};
	double psi(double) const;
	void get_UL_extrapolated (vector<double>&, int) const;
	void get_UR_extrapolated (vector<double>&, int) const;
	void compute_residual(vector<vector<double> >&) const;
};
#endif