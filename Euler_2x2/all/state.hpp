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
public:
	// Constructors
	state()=default;
	state(vector<vector<double> > u, double g) : gamma(g), U(u), CD(false) {};
	state(vector<double>, vector<double>, vector<double>, vector<double>, double, bool=false);
	virtual ~state() = default;
	
	// Set and get members
	inline int get_size() const {return U[0].size();};
	inline int get_dimension() const {return U.size();};
	inline void set_U(const vector<vector<double> >& u) {U=u;return;};
	inline void get_U(vector<vector<double> >& u) const {u=U;return;};
	void get_U(vector<double>& u, int i) const;
	inline void set_gamma(const double g) {gamma = g; return;};
	inline double get_gamma() const {return gamma;};
	inline void set_CD(const bool c) {CD=c;};
	inline bool get_CD() const {return CD;};
	inline void set_sigma(const vector<double>& s) {sigma = s; if(!CD) cerr << "Error: CD set to false" << endl;};
	inline void get_sigma(vector<double>& s) {s = sigma;};
	void print_physical(const string&, ios_base::openmode mode = ios_base::out, int = 15);
	// Methods
	virtual void compute_U_star(const vector<double>&, const vector<double>&, vector<double>&) const = 0;
	virtual void compute_residual(vector<vector<double> >&) const = 0;
	virtual void compute_lambda(vector<double>&) const = 0;
	
};
#endif