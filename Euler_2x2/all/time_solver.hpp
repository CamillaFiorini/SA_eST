#ifndef __TIME_SOLVER_HPP__
#define __TIME_SOLVER_HPP__
#include"state.hpp"
#include"mesh.hpp"
#include"utilities.hpp"
#include<string>

class time_solver
{
protected:
	double start_time;
	double end_time;
	int order;
	mesh m;
	double cfl;
    string path;
public:
	time_solver()=default;
    time_solver(double t, double T, int o, mesh M, double CFL) : start_time(t), end_time(T), order(o), m(M), cfl(CFL), path() {};
    time_solver(double t, double T, int o, mesh M, double CFL, string p) : start_time(t), end_time(T), order(o), m(M), cfl(CFL), path(p) {};
	~time_solver() = default;
	
	inline void set_order(int o) {order = o; return;};
	inline int get_order() const {return order;};
	inline void set_cfl(double c) {cfl = c; return;};
	inline double get_cfl() const {return cfl;};
    inline void set_path(string p) {path = p;};
	
	double solve(state&, bool=false);
};
#endif
