#include"roe_II.hpp"
double roe_II::psi(double r) const
{
	if(r < 1e-15)
		return 0;
	
	//return 1;
	return 0.5*(1 + fabs(r)/r)*min(1., fabs(r));
	//return max(max(0.,min(2.,r)),min(1.,2*r));
	//return (r*r+r)/(1+r*r);
};

void roe_II::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(4);
	int N = this->get_size();
	int left = max(i-1,0);
	int Lleft = max(i-2, 0);
	int right = min(i, N-1);
	
	vector<double> rL(4), irL(4);
	
	for (int k = 0; k < 4; ++k)
	{
		if (fabs(U[k][left]-U[k][Lleft]) < 1e-15 || fabs(U[k][right]-U[k][left]) < 1e-15)
		{
			rL[k] = 0;
			irL[k] = 0;
		}
		else
		{
			rL[k] = (U[k][right]-U[k][left])/(U[k][left]-U[k][Lleft]);
			irL[k] = 1./rL[k];
		}
		UL[k] = U[k][left] +( this->psi(rL[k])*(1-kappa)/4*(U[k][left] - U[k][Lleft]) + this->psi(irL[k])*(1+kappa)/4*(U[k][right] - U[k][left]));
	}/*
	  for (int k = 0; k < 4; ++k)
	  {
	  double a = U[k][right] - U[k][left];
	  double b = U[k][left] - U[k][Lleft];
	  if (a*b < 1e-10 )
			UL[k] = U[k][left];
	  else
			UL[k] = U[k][left] + 0.5*fabs(a)/a*min(fabs(a), fabs(b));
	  }*/
	return;
};

void roe_II::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(4);
	int N = this->get_size();
	int left = max(i-1,0);
	int right = min(i, N-1);
	int Rright = min(i+1, N-1);
	
	vector<double> rR(4), irR(4);
	
	for (int k = 0; k < 4; ++k)
	{
		if (fabs(U[k][Rright]-U[k][right]) < 1e-15 || fabs((U[k][right]-U[k][left])) < 1e-15)
		{
			irR[k] = 0;
			rR[k] = 0;
		}
		else
		{
			rR[k] = (U[k][right]-U[k][left])/(U[k][Rright]-U[k][right]);
			irR[k] = 1./rR[k];
		}
		UR[k] = U[k][right] +1.*(- this->psi(irR[k])*(1+kappa)/4*(U[k][right] - U[k][left]) - this->psi(rR[k])*(1-kappa)/4*(U[k][Rright] - U[k][right]));
	}
	/*
	 for (int k = 0; k < 4; ++k)
	 {
		double a = U[k][Rright] - U[k][right];
		double b = U[k][right] - U[k][left];
		if (a*b < 1e-10 )
	 UR[k] = U[k][right];
		else
	 UR[k] = U[k][right] - 0.5*fabs(a)/a*min(fabs(a), fabs(b));
	 }*/
	return;
};

void roe_II::compute_residual(vector<vector<double> >&) const
{
	return;
};
