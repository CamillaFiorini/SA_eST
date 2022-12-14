#ifndef FLUX_HPP
#define FLUX_HPP

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include"state.hpp"

using namespace std;

class flux
{
public:
	flux() {};
//	void compute(state, vector<vector<double> >&);
	void residual(state, vector<vector<double> >&, vector<vector<double> >&);
	void source(state, vector<vector<double> >&);
	void detector_s1(state, vector<int>&);
	void detector_s2(state, vector<int>&);
};

#endif