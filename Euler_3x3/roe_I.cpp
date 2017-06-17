#include"roe_I.hpp"

void roe_I::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(D);
	int left = max(i-1, 0);
	for (int k = 0; k < D; ++k)
	{
		UL[k] = U[k][left];
	}
	return;
};

void roe_I::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(D);
	int N = U[0].size();
	int right = min(i, N-1);
	for (int k = 0; k < D; ++k)
	{
		UR[k] = U[k][right];
	}
	return;
};