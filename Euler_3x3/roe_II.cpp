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

// returns the Left value of U at the interface i
void roe_II::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(D);
	int N = this->get_size();
	int left = max(i-1,0);
	int Lleft = max(i-2, 0);
	int right = min(i, N-1);
	/*
	vector<double> rL(D), irL(D);
	
	for (int k = 0; k < D; ++k)
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
	for (int k = 0; k < D; ++k)
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

// returns the Right value of U at the interface i
void roe_II::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(D);
	int N = this->get_size();
	int left = max(i-1,0);
	int right = min(i, N-1);
	int Rright = min(i+1, N-1);
	/*
	vector<double> rR(D), irR(D);
	
	for (int k = 0; k < D; ++k)
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
	 for (int k = 0; k < D; ++k)
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
	R.resize(D);
	for (int k = 0; k < D; ++k)
		R[k].assign(N,0);
	vector<double> UL(D), UR(D), F(D/2), UL_star(D), UR_star(D), s_ULstar(D/2), s_URstar(D/2), UL_cell(D), UR_cell(D), s_ULstar_cell(D/2), s_URstar_cell(D/2), UL_star_cell(D),UR_star_cell(D);
	for (int i = 0; i < N+1; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		double lambda1 = this->compute_lambda1(UL, UR);
		double lambda2 = this->compute_lambda2(UL, UR);
		double lambda3 = this->compute_lambda3(UL, UR);
		
		if(CD)
		{
			this->compute_U_star(UL, UR, UL_star, UR_star);
			for (int k = 0; k < D; ++k)
			{
				if(i < N)
				{
					R[k][i] += max(lambda1, 0.0)*UL[k]*(fabs(sigma[i])<1e-10) + (max(lambda2, sigma[i])-max(lambda1, sigma[i]))*UL_star[k] + (max(lambda3,sigma[i])-max(lambda2, sigma[i]))*UR_star[k] - max(lambda3,sigma[i])*UR[k];
				}
				
				if(i > 0)
				{
					R[k][i-1] += min(lambda1, sigma[i])*UL[k] + (min(lambda2, sigma[i]) - min(lambda1, sigma[i]))*UL_star[k] + (min(lambda3,sigma[i])-min(lambda2, sigma[i]))*UR_star[k] - min(lambda3, 0.0)*UR[k]*(fabs(sigma[i])<1e-10);
				}
			}
			if(i > 0)
			{
				UR_cell = UL;
				this->compute_U_star(UL_cell, UR_cell, UL_star_cell, UR_star_cell);
				double lambda1_cell = this->compute_lambda1(UL_cell, UR_cell);
				double lambda2_cell = this->compute_lambda2(UL_cell, UR_cell);
				double lambda3_cell = this->compute_lambda3(UL_cell, UR_cell);
				for (int k = 0; k < D; ++k)
					R[k][i-1] += lambda1_cell*(UL_cell[k] - UL_star_cell[k]) + lambda2_cell*(UL_star_cell[k]-UR_star_cell[k]) + lambda3_cell*(UR_star_cell[k]-UR_cell[k]);
			}
			UL_cell = UR;
		}
		else
		{
			this->compute_U_star(UL, UR, UL_star, UR_star);
			for (int k = 0; k < D; ++k)
			{
				if(i < N)
				{
					R[k][i] += max(lambda1, 0.0)*(UL[k] - UL_star[k]) + max(lambda2, 0.0)*(UL_star[k]-UR_star[k]) + max(lambda3,0.0)*(UR_star[k]-UR[k]);
				}
				
				if(i > 0)
				{
					R[k][i-1] -= min(lambda1, 0.0)*(UL_star[k] - UL[k]) + min(lambda2, 0.0)*(UR_star[k] - UL_star[k]) + min(lambda3,0.0)*(UR[k] - UR_star[k]);
				}
			}
			if(i > 0)
			{
				UR_cell = UL;
				this->compute_U_star(UL_cell, UR_cell, UL_star_cell, UR_star_cell);
				double lambda1_cell = this->compute_lambda1(UL_cell, UR_cell);
				double lambda2_cell = this->compute_lambda2(UL_cell, UR_cell);
				double lambda3_cell = this->compute_lambda3(UL_cell, UR_cell);
				for (int k = 0; k < D; ++k)
					R[k][i-1] += lambda1_cell*(UL_cell[k] - UL_star_cell[k]) + lambda2_cell*(UL_star_cell[k]-UR_star_cell[k]) + lambda3_cell*(UR_star_cell[k]-UR_cell[k]);
			}
			UL_cell = UR;
		}
	}
	return;
};