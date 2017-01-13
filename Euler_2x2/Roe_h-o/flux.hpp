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
	void residual(state, vector<vector<double> >&, vector<vector<double> >&);
	void detector_s1(state, vector<int>&, double=0.0);
	void detector_s2(state, vector<int>&, double=0.0);
	
};

#endif