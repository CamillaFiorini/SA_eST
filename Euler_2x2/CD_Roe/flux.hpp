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
	void source(state, vector<vector<double> >&);
};

#endif