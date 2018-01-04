#ifndef ROE_I_HPP
#define ROE_I_HPP

#include<vector>
#include"roe.hpp"

using namespace std;

class roe_I:public roe
{
public:
	// Constructors
	roe_I(const vector<vector<double> >& u, double g) : roe(u,g) {};
	roe_I(const vector<double>&a, const vector<double>& b, const vector<double>& c, const vector<double>& d, const vector<double>& e, const vector<double>& f, double g,const vector<double>& h, const vector<double>& dh) : roe(a,b,c,d,e,f,g,h,dh) {};
	roe_I(const vector<double>&a, const vector<double>& b, const vector<double>& c, const vector<double>& d, const vector<double>& e, const vector<double>& f, double g,const vector<double>& h, const vector<double>& dh, const vector<double>& sh, const vector<double>& dsh) : roe(a,b,c,d,e,f,g,h,dh,sh,dsh) {};
	// Methods
	void get_UL_extrapolated (vector<double>&, int) const;
	void get_UR_extrapolated (vector<double>&, int) const;
};
#endif