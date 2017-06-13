#include"roe_I.hpp"

void roe_I::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(4);
	int left = max(i-1, 0);
	for (int k = 0; k < 4; ++k)
	{
		UL[k] = U[k][left];
	}
	return;
};

void roe_I::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(4);
	int N = U[0].size();
	int right = min(i, N-1);
	for (int k = 0; k < 4; ++k)
	{
		UR[k] = U[k][right];
	}
	return;
};