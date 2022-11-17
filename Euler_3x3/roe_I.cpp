#include"roe_I.hpp"

void roe_I::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(D);
	if(i > 0)
		for (int k = 0; k < D; ++k)
			UL[k] = U[k][i-1];
	if(i == 0)
		for (int k = 0; k < D; ++k)
			UL[k] = U[k][i];
	
	if(i == 0 && (bc_L[0] || bc_L[1] || bc_L[2]))
	{
		vector<double> Ui, W; //conservative and physical variables interior cell Ui = (rho, rho u, rho E), W = (rho, u, p)
		this->get_W(W,0);
		this->get_U(Ui,0);
		double H, p_tot, p, s_H, s_p_tot, s_p;
		
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
		
		UL[0] = (p/(gamma-1)+p_tot)/H;
		UL[3] = (s_p/(gamma-1)+s_p_tot - UL[0]*s_H)/H;
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
		UL[2] = p/(gamma-1)+0.5*UL[1]*UL[1]/UL[0]; // rhoE = p/(gamma-1)+ 0.5rho*u^2
		UL[5] = s_p/(gamma-1)+UL[1]*UL[4]/UL[0]-0.5*UL[3]*UL[1]*UL[1]/UL[0]/UL[0];
		
	}
	return;
};

void roe_I::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(D);
	int N = U[0].size();
	double H, p_tot, p, s_H, s_p_tot, s_p;

	if(i < N)
		for (int k = 0; k < D; ++k)
			UR[k] = U[k][i];
	if(i == N)
		for (int k = 0; k < D; ++k)
			UR[k] = U[k][i-1];
	
	if(i == N && (bc_R[0] || bc_R[1] || bc_R[2]))
	{
		vector<double> Ui, W; //conservative and physical variables interior cell
		this->get_W(W,N-1);
		this->get_U(Ui,N-1);
		
		if(bc_R[0])
		{
			H = VR_inf[0];
			s_H = VR_inf[3];
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
			s_p =  W[5];
		}
		
		UR[0] = (p/(gamma-1)+p_tot)/H;
		UR[3] = (s_p/(gamma-1)+s_p_tot-UR[0]*s_H)/H;
		if (p_tot >= p)
		{
			UR[1] = sqrt(2*UR[0]*(p_tot-p));
			UR[4] = 1/UR[1]*(UR[3]*(p_tot-p) + UR[0]*(s_p_tot-s_p));
		}
		else
		{
			UR[1] = 0;
			UR[4] = 0;
		}
		UR[2] = p/(gamma-1)+0.5*UR[1]*UR[1]/UR[0];
		UR[5] = s_p/(gamma-1)+UR[1]*UR[4]/UR[0]-0.5*UR[3]*UR[1]*UR[1]/UR[0]/UR[0];
	}
	return;
};
