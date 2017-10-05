#include"roe_I.hpp"

void roe_I::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(D);
	if(i > 0)
	{
		int left = max(i-1,0);
		for (int k = 0; k < D; ++k)
		{
			UL[k] = U[k][left];
		}
	}
	else
	{
		vector<double> Ui, W; //physical variables interior cell = (rho, u, p)
		this->get_W(W,0);
		this->get_U(Ui,0);
		double H, p_tot, p;
		
		if(bc_L[0])
			H = VL_inf[0];
		else
			H = this->compute_H(Ui);
		
		if(bc_L[1])
			p_tot = VL_inf[1];
		else
			p_tot = W[2]+0.5*W[0]*W[1]*W[1];
		
		if (bc_L[2])
			p = VL_inf[2];
		else
			p = W[2];
		
		UL[0] = (p/(gamma-1)+p_tot)/H;
		if (p_tot > p)
			UL[1] = sqrt(2*UL[0]*(p_tot-p));
		else
			UL[1] = 0;
		UL[2] = p/(gamma-1)+p_tot-p;
	}
	return;
};

void roe_I::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(D);
	int N = U[0].size();
	
	if(i < N)
	{
		int right = min(i, N-1);
		for (int k = 0; k < D; ++k)
		{
			UR[k] = U[k][right];
		}
	}
	else
	{
		vector<double> Ui, W; //conservative and physical variables interior cell
		this->get_W(W,N-1);
		this->get_U(Ui,N-1);
		double H, p_tot, p;
		
		if(bc_R[0])
			H = VR_inf[0];
		else
			H = this->compute_H(Ui);
		
		if(bc_R[1])
			p_tot = VR_inf[1];
		else
			p_tot = W[2]+0.5*W[0]*W[1]*W[1];
		
		if (bc_R[2])
			p = VR_inf[2];
		else
			p = W[2];
		
		UR[0] = (p/(gamma-1)+p_tot)/H;
		UR[1] = sqrt(2*UR[0]*(p_tot-p));
		UR[2] = p/(gamma-1)+p_tot-p;
	}
	return;
};