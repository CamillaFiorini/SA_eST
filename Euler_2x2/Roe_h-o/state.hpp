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
	
	// Functions that return the extrapolated values at the interface i
	void get_UL_extrapolated (vector<double>&, int, double=1./3.);
	void get_UR_extrapolated (vector<double>&, int, double=1./3.);

	double psi(double);
	
	void compute_lambdaR(vector<double>&);
	void compute_s_lambdaR(vector<double>&);
	
//	void detector_s1(vector<int>&, double=0.0);
//	void detector_s2(vector<int>&, double=0.0);
};

#endif