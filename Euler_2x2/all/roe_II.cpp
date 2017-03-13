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
	/*
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
	}*/
	for (int k = 0; k < 4; ++k)
	{
		double a = U[k][right] - U[k][left];
		double b = U[k][left] - U[k][Lleft];
		
		if (a*b < 1e-10 )
			UL[k] = U[k][left];
		else
			UL[k] = U[k][left] + 0.5*fabs(a)/a*min(fabs(a), fabs(b));
	}
	return;
};

void roe_II::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(4);
	int N = this->get_size();
	int left = max(i-1,0);
	int right = min(i, N-1);
	int Rright = min(i+1, N-1);
	/*
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
	*/
	 for (int k = 0; k < 4; ++k)
	 {
		double a = U[k][Rright] - U[k][right];
		double b = U[k][right] - U[k][left];
		if (a*b < 1e-10 )
			UR[k] = U[k][right];
		else
			UR[k] = U[k][right] - 0.5*fabs(a)/a*min(fabs(a), fabs(b));
	 }
	return;
};

void roe_II::compute_residual(vector<vector<double> >& R) const
{
	int N = this->get_size();
	R.resize(4);
	for (int k = 0; k < 4; ++k)
		R[k].assign(N,0);
	
	vector<double> UL(4), UR(4), Ustar(4), UL_cell(4), UR_cell(4), Ustar_cell(4);
	double lambda, lambda_cell;
	
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		
		lambda = this->compute_lambda(UL, UR);
		this->compute_U_star(UL, UR, Ustar);
		
		for (int k=0; k<4; ++k)
		{
			if (i < N)
				R[k][i] += lambda*(Ustar[k]-UR[k]);
			
			if (i > 0)
				R[k][i-1] += lambda*(Ustar[k]-UL[k]);
		}
		
		if(i > 0)
		{
			UR_cell = UL;
			this->compute_U_star(UL_cell, UR_cell, Ustar_cell);
			lambda_cell = this->compute_lambda(UL_cell, UR_cell);
			for (int k = 0; k < 4; ++k)
				R[k][i-1] += lambda_cell*(2*Ustar_cell[k] - UL_cell[k] - UR_cell[k]);
		}
		UL_cell = UR;
		
	}

	
	return;
};
