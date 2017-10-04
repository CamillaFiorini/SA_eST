#ifndef ROE_I_HPP
#define ROE_I_HPP

#include<vector>
#include"roe.hpp"

using namespace std;

class roe_I:public roe
{
public:
	// Constructors
	roe_I(vector<vector<double> > u, double g) : roe(u,g) {};
	roe_I(vector<double>a, vector<double>b, vector<double>c, vector<double>d, vector<double>e, vector<double>f, double g, vector<double>h, vector<double> dh) : roe(a,b,c,d,e,f,g,h,dh) {};
	// Methods
	void get_UL_extrapolated (vector<double>&, int) const;
	void get_UR_extrapolated (vector<double>&, int) const;
};
#endif