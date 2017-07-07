#ifndef ROE_HPP
#define ROE_HPP

#include<vector>
#include<algorithm>
#include"state.hpp"

using namespace std;

class roe:public state
{
private:
	bool sens_hllc;
public:
	// Constructors
	roe(vector<vector<double> > u, double g) : state(u,g), sens_hllc(false) {};
	roe(vector<double>a, vector<double>b, vector<double>c, vector<double>d, vector<double>e, vector<double>f, double g) : state(a,b,c,d,e,f,g), sens_hllc(false) {};
	// Methods
	inline void set_sens_hllc(bool c) {sens_hllc = c;};
	inline bool get_sens_hllc() {return sens_hllc;};
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
	int detector_s3(const vector<double>&, const vector<double>&, double=1e-10) const;
	int detector_c(const vector<double>&, const vector<double>&, double=1e-10) const;
	void detector_s1(vector<int>&, double=1e-10) const;
	void detector_s3(vector<int>&, double=1e-10) const;
	void detector_c(vector<int>&, double=1e-10) const;
	void compute_flux(const vector<double>&, const vector<double>&, vector<double>&, vector<double>&, vector<double>&, int i = 0) const;
	void compute_U_star(const vector<double>&, const vector<double>&, vector<double>&, vector<double>&) const;
	virtual void get_UL_extrapolated (vector<double>&, int) const = 0;
	virtual void get_UR_extrapolated (vector<double>&, int) const = 0;
	virtual void compute_residual(vector<vector<double> >&) const;
};
#endif