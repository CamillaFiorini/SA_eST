#ifndef STATE_HPP
#define STATE_HPP

#include<vector>
#include<cmath>
#include<iostream>

using namespace std;

class state
{
protected:
	double gamma;
	vector<vector<double> > U;
public:
	// Constructors
	state()=default;
	state(vector<vector<double> > u, double g) : gamma(g), U(u) {};
	state(vector<double>, vector<double>, vector<double>, vector<double>, double);
	// ToutDoux! Destructor
	//virtual ~state() = default;
	
	// Set and get members
	inline int get_size() const {return U[0].size();};
	inline void set_U(const vector<vector<double> >& u) {U=u;return;};
	inline void get_U(vector<vector<double> >& u) const {u=U;return;};
	inline void set_gamma(const double g) {gamma = g; return;};
	inline double get_gamma() const {return gamma;};
	
	// Methods
	virtual void compute_lambda(vector<double>&) const = 0;
	virtual void compute_U_star(const vector<double>&, const vector<double>&, vector<double>&) const = 0;
	virtual void compute_residual(vector<vector<double> >&) const = 0;
	
};
#endif