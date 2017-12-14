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
	vector<double> Ui, W; //conservative and physical variables interior cell Ui = (rho, rho u, rho E), W = (rho, u, p)
	double H, p_tot, p, s_H, s_p_tot, s_p;
	if(i < 2) // for the first two interfaces I compute values "outside" the domain (depending on the BC)
	{
		this->get_W(W,0);
		this->get_U(Ui,0);
		if(bc_L[0])
		{
			H = VL_inf[0];
			s_H = VL_inf[3];
		}
		else
		{
			H = this->compute_H(Ui);
			s_H = this->compute_s_H(Ui);
		}
		
		if(bc_L[1])
		{
			p_tot = VL_inf[1];
			s_p_tot = VL_inf[4];
		}
		else
		{
			p_tot = W[2]+0.5*W[0]*W[1]*W[1];
			s_p_tot = W[5]+0.5*W[3]*W[1]*W[1]+W[0]*W[1]*W[4];
		}
		
		if (bc_L[2])
		{
			p = VL_inf[2];
			s_p = VL_inf[5];
		}
		else
		{
			p = W[2];
			s_p = W[5];
		}
	}
	if(i == 0) // at the first interface I impose the value outside the domain
	{
		UL[0] = (p/(gamma-1)+p_tot)/H;
		if (p_tot > p)
		{
			UL[1] = sqrt(2*UL[0]*(p_tot-p));
			UL[4] = 1/UL[1]*(UL[3]*(p_tot-p) + UL[0]*(s_p_tot-s_p));
		}
		else
		{
			UL[1] = 0;
			UL[4] = 0;
		}
		UL[2] = p/(gamma-1)+p_tot-p;
		UL[3] = (s_p/(gamma-1)+s_p_tot)/H - UL[0]*s_H/H;
		UL[5] = s_p/(gamma-1)+s_p_tot-s_p;
	}
	if(i == 1) // at the second interface I use the value outside the domain as a Lleft value
	{
		vector<double> ULleft(D);
		ULleft[0] = (p/(gamma-1)+p_tot)/H;
		if (p_tot > p)
		{
			ULleft[1] = sqrt(2*ULleft[0]*(p_tot-p));
			ULleft[4] = 1/ULleft[1]*(ULleft[3]*(p_tot-p) + ULleft[0]*(s_p_tot-s_p));
		}
		else
		{
			ULleft[1] = 0;
			ULleft[4] = 0;
		}
		ULleft[2] = p/(gamma-1)+p_tot-p;
		ULleft[3] = (s_p/(gamma-1)+s_p_tot)/H - ULleft[0]*s_H/H;
		ULleft[5] = s_p/(gamma-1)+s_p_tot-s_p;
		
		for (int k = 0; k < D; ++k)
		{
			double a = U[k][right] - U[k][left];
			double b = U[k][left] - ULleft[k];
			
			if (a*b < 1e-10 )
				UL[k] = U[k][left];
			else
				UL[k] = U[k][left] + 0.5*fabs(a)/a*min(fabs(a), fabs(b));
		}
	}
	if(i > 1) // internal interfaces
	{
		for (int k = 0; k < D; ++k)
		{
			double a = U[k][right] - U[k][left];
			double b = U[k][left] - U[k][Lleft];
			
			if (a*b < 1e-10 )
				UL[k] = U[k][left];
			else
				UL[k] = U[k][left] + 0.5*fabs(a)/a*min(fabs(a), fabs(b));
		}
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
	vector<double> Ui, W; //conservative and physical variables interior cell
	double H, p_tot, p, s_H, s_p_tot, s_p;
	if(i > N-2) // for the last two interfaces I compute the "outside" values (according to the BC)
	{
		this->get_W(W,N-1);
		this->get_U(Ui,N-1);
		if(bc_R[0])
		{
			H = VR_inf[0];
			s_H = VR_inf[0];
		}
		else
		{
			H = this->compute_H(Ui);
			s_H = this->compute_s_H(Ui);
		}
		if(bc_R[1])
		{
			p_tot = VR_inf[1];
			s_p_tot = VR_inf[4];
		}
		else
		{
			p_tot = W[2]+0.5*W[0]*W[1]*W[1];
			s_p_tot = W[5]+0.5*W[3]*W[1]*W[1]+W[0]*W[1]*W[4];
		}
		
		if (bc_R[2])
		{
			p = VR_inf[2];
			s_p = VR_inf[5];
		}
		else
		{
			p = W[2];
			s_p = W[5];
		}
	}
	if (i == N) // last interface: the outside values are directly imposed
	{
		UR[0] = (p/(gamma-1)+p_tot)/H;
		UR[1] = sqrt(2*UR[0]*(p_tot-p));
		UR[2] = p/(gamma-1)+p_tot-p;
		UR[3] = (s_p/(gamma-1)+s_p_tot)/H - UR[0]*s_H/H;;
		UR[4] = 1/UR[1]*(UR[3]*(p_tot-p) + UR[0]*(s_p_tot-s_p));
		UR[5] = s_p/(gamma-1)+s_p_tot-s_p;
	}
	if(i == N-1) // second to last interface: outside values used as Rright values
	{
		vector<double> URright(D);
		URright[0] = (p/(gamma-1)+p_tot)/H;
		URright[1] = sqrt(2*URright[0]*(p_tot-p));
		URright[2] = p/(gamma-1)+p_tot-p;
		URright[3] = (s_p/(gamma-1)+s_p_tot)/H - URright[0]*s_H/H;;
		URright[4] = 1/URright[1]*(URright[3]*(p_tot-p) + URright[0]*(s_p_tot-s_p));
		URright[5] = s_p/(gamma-1)+s_p_tot-s_p;
		
		for (int k = 0; k < D; ++k)
		{
			double a = URright[k] - U[k][right];
			double b = U[k][right] - U[k][left];
			if (a*b < 1e-10 )
				UR[k] = U[k][right];
			else
				UR[k] = U[k][right] - 0.5*fabs(a)/a*min(fabs(a), fabs(b));
		}
	}
	if (i < N-1) // internal interfaces
	{
		for (int k = 0; k < D; ++k)
		{
			double a = U[k][Rright] - U[k][right];
			double b = U[k][right] - U[k][left];
			if (a*b < 1e-10 )
				UR[k] = U[k][right];
			else
				UR[k] = U[k][right] - 0.5*fabs(a)/a*min(fabs(a), fabs(b));
		}
	}
	return;
};

void roe_II::compute_residual(vector<vector<double> >& R) const
{
	int N = this->get_size();
	R.resize(D);
	for (int k = 0; k < D; ++k)
		R[k].assign(N,0);
	vector<double> UL(D), UR(D), F(D/2), UL_star(D), UR_star(D), U_SL(D), U_SR(D), s_ULstar(D/2), s_URstar(D/2), UL_cell(D), UR_cell(D), s_ULstar_cell(D/2), s_URstar_cell(D/2), UL_star_cell(D),UR_star_cell(D), lambda(3);
	for (int i = 0; i < N+1; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		this->compute_U_star(UL, UR, UL_star, UR_star, lambda, U_SL, U_SR);
		if(CD)
		{
			for (int k = 0; k < D; ++k)
			{
				if(i < N)
				{
					R[k][i] += max(lambda[0], 0.0)*UL[k]*(fabs(sigma[i])<1e-10) + (max(lambda[1], sigma[i])-max(lambda[0], sigma[i]))*UL_star[k] + (max(lambda[2],sigma[i])-max(lambda[1], sigma[i]))*UR_star[k] - max(lambda[2],sigma[i])*UR[k];
				}
				
				if(i > 0)
				{
					R[k][i-1] += min(lambda[0], sigma[i])*UL[k] + (min(lambda[1], sigma[i]) - min(lambda[0], sigma[i]))*UL_star[k] + (min(lambda[2],sigma[i])-min(lambda[1], sigma[i]))*UR_star[k] - min(lambda[2], 0.0)*UR[k]*(fabs(sigma[i])<1e-10);
				}
			}
			if(i > 0)
			{
				UR_cell = UL;
				vector<double> lambda_cell, U_SL_cell, U_SR_cell;
				this->compute_U_star(UL_cell, UR_cell, UL_star_cell, UR_star_cell, lambda_cell, U_SL_cell, U_SR_cell);
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
			for (int k = 0; k < D; ++k)
			{
				if(i < N)
				{
					R[k][i] += max(lambda[0], 0.0)*(UL[k] - UL_star[k]) + max(lambda[1], 0.0)*(UL_star[k]-UR_star[k]) + max(lambda[2],0.0)*(UR_star[k]-UR[k]);
				}
				
				if(i > 0)
				{
					R[k][i-1] -= min(lambda[0], 0.0)*(UL_star[k] - UL[k]) + min(lambda[1], 0.0)*(UR_star[k] - UL_star[k]) + min(lambda[2],0.0)*(UR[k] - UR_star[k]);
				}
			}
			if(i > 0)
			{
				UR_cell = UL;
				vector<double> lambda_cell(3), U_SL_cell, U_SR_cell;
				this->compute_U_star(UL_cell, UR_cell, UL_star_cell, UR_star_cell, lambda_cell, U_SL_cell, U_SR_cell);
				double lambda1_cell = this->compute_lambda1(UL_cell, UR_cell);
				double lambda2_cell = this->compute_lambda2(UL_cell, UR_cell);
				double lambda3_cell = this->compute_lambda3(UL_cell, UR_cell);
				for (int k = 0; k < D; ++k)
					R[k][i-1] += lambda1_cell*(UL_cell[k] - UL_star_cell[k]) + lambda2_cell*(UL_star_cell[k]-UR_star_cell[k]) + lambda3_cell*(UR_star_cell[k]-UR_cell[k]);
			}
			UL_cell = UR;
		}
		if(i < N)
		{
			// adding (dx(h))(P-F)
			vector<double> Wi(D), Ui(D), Fi(D);
			this->get_U(Ui,i);
			this->get_W(Wi,i);
			double p (Wi[2]), s_p(Wi[5]);
			this->flux(Ui,Fi);
			R[0][i] += -Fi[0]*delta_h[i]/h[i];
			R[1][i] += (p - Fi[1])*delta_h[i]/h[i];
			R[2][i] += -Fi[2]*delta_h[i]/h[i];
			R[3][i] += -Fi[3]*delta_h[i]/h[i];
			R[4][i] += (s_p - Fi[4])*delta_h[i]/h[i];
			R[5][i] += -Fi[5]*delta_h[i]/h[i];
		}
	}
	return;
};