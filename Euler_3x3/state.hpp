#ifndef STATE_HPP
#define STATE_HPP

#include<vector>
#include<cmath>
#include<iostream>
#include<string>
#include<fstream>
#include<omp.h>

using namespace std;

class state
{
protected:
	double gamma;
	vector<vector<double> > U;
	bool CD;
	vector<double> sigma;
	int D;
public:
	// Constructors
	state()=default;
	state(const vector<vector<double> >& u, double g) : gamma(g), U(u), CD(false), D(u.size()) {};
	state(const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, double, bool=false);
	virtual ~state() = default;
	
	// Set and get members
	inline int get_size() const {return U[0].size();};
	inline int get_dimension() const {return U.size();};
	inline void set_U(const vector<vector<double> >& u) {U=u; D=u.size(); return;};
	inline void get_U(vector<vector<double> >& u) const {u=U;return;};
	void get_U(vector<double>& u, int i) const;
	void get_W(vector<vector<double> >& a) const;
	void flux(const vector<double>&, vector<double>&) const;
	inline void set_gamma(const double g) {gamma = g; return;};
	inline double get_gamma() const {return gamma;};
	inline void set_CD(const bool c) {CD=c;};
	inline bool get_CD() const {return CD;};
	inline void set_sigma(const vector<double>& s) {sigma = s; if(!CD) cerr << "Error: CD set to false" << endl;};
	inline void get_sigma(vector<double>& s) {s = sigma;};
	void print_conservative(const string&, ios_base::openmode mode = ios_base::out, int = 15);
	void print_physical(const string&, ios_base::openmode mode = ios_base::out, int = 15);
	
	// Virtual methods
	virtual void compute_residual(vector<vector<double> >&) const = 0;
	virtual double compute_maxvel() const = 0;
	virtual double compute_lambda1(const vector<double>&, const vector<double>&) const = 0;
	virtual void compute_lambda1(vector<double>&) const = 0;
	virtual double compute_lambda2(const vector<double>&, const vector<double>&) const = 0;
	virtual void compute_lambda2(vector<double>&) const = 0;
	virtual double compute_lambda3(const vector<double>&, const vector<double>&) const = 0;
	virtual void compute_lambda3(vector<double>&) const = 0;
	virtual void detector_s1(vector<int>&, double=1e-10) const = 0;
	virtual void detector_s3(vector<int>&, double=1e-10) const = 0;
	virtual void detector_c(vector<int>&, double=1e-10) const = 0;
	
};
#endif