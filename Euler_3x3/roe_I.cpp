#include"roe_I.hpp"

void roe_I::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(D);
	if(!bc_L || i > 0)
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
		double H = VL_inf[0]; //this->compute_H(Ui);
		double p_tot = VL_inf[1]; //W[2]+0.5*W[0]*W[1]*W[1];
		double p = W[2]; //VL_inf[2];
		UL[0] = (p/(gamma-1)+p_tot)/H;
		UL[1] = sqrt(2*UL[0]*(p_tot-p));
		UL[2] = p/(gamma-1)+p_tot-p;
	}
	return;
};

void roe_I::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(D);
	int N = U[0].size();
	
	if(!bc_R || i < N)
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
		double H = this->compute_H(Ui); //VR_inf[0];
		double p_tot = W[2]+0.5*W[0]*W[1]*W[1]; //VR_inf[1];
		double p =  VR_inf[2]; // W[2];
		UR[0] = (p/(gamma-1)+p_tot)/H;
		UR[1] = sqrt(2*UR[0]*(p_tot-p));
		UR[2] = p/(gamma-1)+p_tot-p;
	}
	return;
};