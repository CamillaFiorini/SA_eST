#include"roe.hpp"

void roe::compute_lambda1(vector<double>& l) const
{
	int N = U[0].size();
	l.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		l[i] = this->compute_lambda1(UL, UR);
	}
	return;
};

double roe::compute_lambda1(const vector<double>& UL, const vector<double>& UR) const
{
	double utilde = this->compute_utilde(UL, UR);
	double atilde = this->compute_atilde(UL, UR);
	return utilde - atilde;
};

void roe::compute_lambda2(vector<double>& l) const
{
	int N = U[0].size();
	l.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		l[i] = this->compute_lambda2(UL, UR);
	}
	return;
};

double roe::compute_lambda2(const vector<double>& UL, const vector<double>& UR) const
{
	return this->compute_utilde(UL, UR);
};

void roe::compute_lambda3(vector<double>& l) const
{
	int N = U[0].size();
	l.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		l[i] = this->compute_lambda3(UL, UR);
	}
	return;
};

double roe::compute_lambda3(const vector<double>& UL, const vector<double>& UR) const
{
	double utilde = this->compute_utilde(UL, UR);
	double atilde = this->compute_atilde(UL, UR);
	return utilde + atilde;
};

double roe::compute_maxvel() const
{
	vector<double> l1, l2, l3;
	this->compute_lambda1(l1);
	this->compute_lambda2(l2);
	this->compute_lambda3(l3);
	double M = max(fabs(l1[0]), max(fabs(l2[0]), fabs(l3[0])));
	for (unsigned int i = 1; i < l1.size(); ++i)
	{
		if (M < max(fabs(l1[i]), max(fabs(l2[i]), fabs(l3[i]))))
			M = max(fabs(l1[i]), max(fabs(l2[i]), fabs(l3[i])));
	}
	return M;
};

double roe::compute_utilde(const vector<double>& UL, const vector<double>& UR) const
{
	return (UL[1]/sqrt(UL[0]) + UR[1]/sqrt(UR[0]))/( sqrt(UL[0]) + sqrt(UR[0]) );
};

double roe::compute_s_utilde(const vector<double>& UL, const vector<double>& UR) const
{
	double uR = UR[1]/UR[0];
	double uL = UL[1]/UL[0];
	double s_uR = UR[4]/UR[0] - UR[3]*UR[1]/(UR[0]*UR[0]);
	double s_uL = UL[4]/UL[0] - UL[3]*UL[1]/(UL[0]*UL[0]);
	double u_tilde = this->compute_utilde(UL, UR);
	return 0.5*UR[3]*uR/sqrt(UR[0])/(sqrt(UR[0])+sqrt(UL[0])) + 0.5*UL[3]*uL/sqrt(UL[0])/(sqrt(UR[0])+sqrt(UL[0])) + (sqrt(UR[0])*s_uR + sqrt(UL[0])*s_uL)/(sqrt(UR[0])+sqrt(UL[0])) - u_tilde/(sqrt(UR[0])+sqrt(UL[0]))*(0.5*UR[3]/sqrt(UR[0]) + 0.5*UL[3]/sqrt(UL[0]));
};

double roe::compute_Htilde(const vector<double>& UL, const vector<double>& UR) const
{
	double HL = this->compute_H(UL);
	double HR = this->compute_H(UR);
	return (sqrt(UL[0])*HL + sqrt(UR[0])*HR)/(sqrt(UL[0]) + sqrt(UR[0]));
};

double roe::compute_s_Htilde(const vector<double>& UL, const vector<double>& UR) const
{
	double HR = this->compute_H(UR);
	double HL = this->compute_H(UL);
	double s_HR = this->compute_s_H(UR);
	double s_HL = this->compute_s_H(UL);
	double H_tilde = this->compute_Htilde(UL, UR);
	return 0.5*UR[3]*HR/sqrt(UR[0])/(sqrt(UR[0])+sqrt(UL[0])) + 0.5*UL[3]*HL/sqrt(UL[0])/(sqrt(UR[0])+sqrt(UL[0])) + (sqrt(UR[0])*s_HR + sqrt(UL[0])*s_HL)/(sqrt(UR[0])+sqrt(UL[0])) - H_tilde/(sqrt(UR[0])+sqrt(UL[0]))*(0.5*UR[3]/sqrt(UR[0]) + 0.5*UL[3]/sqrt(UL[0]));
};

double roe::compute_atilde(const vector<double>& UL, const vector<double>& UR) const
{
	double utilde = this->compute_utilde(UL, UR);
	double Htilde = this->compute_Htilde(UL, UR);
	return sqrt((gamma-1)*(Htilde - 0.5*utilde*utilde));
};

double roe::compute_s_atilde(const vector<double>& UL, const vector<double>& UR) const
{
	double s_Htilde = this->compute_s_Htilde(UL, UR);
	double u_tilde = this->compute_utilde(UL, UR);
	double s_u_tilde = this->compute_s_utilde(UL, UR);
	double H_tilde = this->compute_Htilde(UL, UR);
	return 0.5*(gamma-1)*(s_Htilde - u_tilde*s_u_tilde)/sqrt((gamma-1)*(H_tilde - 0.5*u_tilde*u_tilde));
};

void roe::compute_alpha_tilde(const vector<double>& UL, const vector<double>& UR, vector<double>& alpha_tilde) const
{
	alpha_tilde.resize(D);
	double utilde = this->compute_utilde(UL, UR);
	double atilde = this->compute_atilde(UL, UR);
	double Htilde = this->compute_Htilde(UL, UR);
	double s_utilde = this->compute_s_utilde(UL, UR);
	double s_atilde = this->compute_s_atilde(UL, UR);
	double s_Htilde = this->compute_s_Htilde(UL, UR);

	alpha_tilde[1] = (gamma-1)/(atilde*atilde)*( (UR[0] - UL[0])*(Htilde - utilde*utilde) + utilde*(UR[1] - UL[1]) - (UR[2] - UL[2]));
	alpha_tilde[0] = 0.5/atilde*( (UR[0] - UL[0])*(utilde+atilde) - (UR[1] - UL[1]) - atilde*alpha_tilde[1] );
	alpha_tilde[2] = (UR[0] - UL[0]) - (alpha_tilde[0] + alpha_tilde[1]);
	
	alpha_tilde[4] =	-2*s_atilde*(gamma-1)/(atilde*atilde*atilde)*( (UR[0] - UL[0])*(Htilde - utilde*utilde) + utilde*(UR[1] - UL[1]) - (UR[2] - UL[2]))
						+(gamma-1)/(atilde*atilde)*( (UR[3] - UL[3])*(Htilde - utilde*utilde) + (UR[0] - UL[0])*(s_Htilde - 2*s_utilde*utilde) + s_utilde*(UR[1] - UL[1]) + utilde*(UR[4] - UL[4]) - (UR[5] - UL[5]));
	
	alpha_tilde[3] =	-s_atilde*0.5/(atilde*atilde)*( (UR[0] - UL[0])*(utilde+atilde) - (UR[1] - UL[1]) - atilde*alpha_tilde[1] )
						+ 0.5/atilde*( (UR[3] - UL[3])*(utilde+atilde) + (UR[0] - UL[0])*(s_utilde+s_atilde) - (UR[4] - UL[4]) - s_atilde*alpha_tilde[1] - atilde*alpha_tilde[4] );
	
	alpha_tilde[5] = (UR[3] - UL[3]) - (alpha_tilde[3] + alpha_tilde[4]);
	return;
};

int roe::detector_s1(const vector<double>& UL, const vector<double>& UR, double threshold) const
{
	if ((UR[0] - UL[0]) > threshold && (UL[1]/UL[0] - UR[1]/UR[0]) > threshold )
		return 1;
	return 0;
};

int roe::detector_s3(const vector<double>& UL, const vector<double>& UR, double threshold) const
{
	if ((UL[0] - UR[0]) > threshold && (UL[1]/UL[0] - UR[1]/UR[0]) > threshold )
		return 1;
	return 0;
};

int roe::detector_c(const vector<double>& UL, const vector<double>& UR, double threshold) const
{
	if ((UL[0] - UR[0]) > threshold && fabs(UL[1]/UL[0] - UR[1]/UR[0]) < threshold )
		return 1;
	return 0;
};

void roe::detector_s1(vector<int>& d1, double threshold) const
{
	int N = U[0].size();
	d1.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		d1[i] = this->detector_s1(UL, UR, threshold);
	}
	return;
};

void roe::detector_s3(vector<int>& d3, double threshold) const
{
	int N = U[0].size();
	d3.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		d3[i] = this->detector_s3(UL, UR, threshold);
	}
	return;
};

void roe::detector_c(vector<int>& c, double threshold) const
{
	int N = U[0].size();
	c.resize(N+1);
	vector<double> UL, UR;
	for (int i = 0; i <= N; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		c[i] = this->detector_c(UL, UR, threshold);
	}
	return;
};

int roe::compute_U_star(const vector<double>& UL, const vector<double>& UR, vector<double>& UL_star, vector<double>& UR_star, vector<double>& lambda, vector<double>& U_SL, vector<double>& U_SR, int i) const
{
	int transonic_raref = 0;
	UL_star.resize(D);
	UR_star.resize(D);
	U_SL.resize(D);
	U_SR.resize(D);
		
	vector<double> r1(D), r3(D), alpha_tilde(D), S(D/2);
	int d1(this->detector_s1(UL, UR, 4e-3)), d3(this->detector_s3(UL, UR, 4e-3));
	this->compute_alpha_tilde(UL, UR, alpha_tilde);
	
	double utilde = this->compute_utilde(UL, UR);
	double atilde = this->compute_atilde(UL, UR);
	double Htilde = this->compute_Htilde(UL, UR);
	double s_atilde = this->compute_s_atilde(UL, UR);
	double s_utilde = this->compute_s_utilde(UL, UR);
	double s_Htilde = this->compute_s_Htilde(UL, UR);
	vector<double> FL(D), FR(D);
	this->flux(UL, FL);
	this->flux(UR, FR);
	lambda.resize(3);
	lambda[0] = utilde-atilde;
	lambda[1] = utilde;
	lambda[2] = utilde+atilde;
	
	r1[0] = 1;
	r1[1] = utilde - atilde;
	r1[2] = Htilde - utilde*atilde;
	r1[3] = 0;
	r1[4] = s_utilde - s_atilde;
	r1[5] = s_Htilde - s_utilde*atilde - utilde*s_atilde;
	
	r3[0] = 1;
	r3[1] = utilde + atilde;
	r3[2] = Htilde + utilde*atilde;
	r3[3] = 0;
	r3[4] = s_utilde + s_atilde;
	r3[5] = s_Htilde + s_utilde*atilde + utilde*s_atilde;

	for (int k = 0; k < D/2; ++k)
	{
		UL_star[k] = UL[k] + alpha_tilde[0]*r1[k];
		UR_star[k] = UR[k] - alpha_tilde[2]*r3[k];
	}
	// Checking for left transonic rarefactions
	vector<double> WL, WL_star;
	this->conservative2physical(UL, WL);
	this->conservative2physical(UL_star, WL_star);
	double lambda1L = WL[1] - sqrt(gamma*WL[2]/WL[0]);
	double lambda1L_star = WL_star[1] - sqrt(gamma*WL_star[2]/WL_star[0]);
	//left transonic rarefaction: entropy fix
	if(lambda1L < 0 && lambda1L_star > 0)
	{
		transonic_raref = 1;
		for (int k = 0; k < D; ++k)
		{
			U_SL[k] = ((lambda[0]-lambda1L)*UL[k] + (lambda1L_star-lambda[0])*UL_star[k])/(lambda1L_star-lambda1L);
		}
	}
	for (int k = 0; k < D/2; ++k)
	{
		S[k] = (s_utilde - s_atilde)*(UL_star[k] - UL[k])*d1 + s_utilde*(UR_star[k] - UL_star[k]) + (s_utilde + s_atilde)*(UR[k] - UR_star[k])*d3;
		if (!sens_hllc)
		{
			UR_star[k+3] = 1.0/(lambda[2]-lambda[0])*(lambda[2]*UR[k+3] - lambda[0]*UL[k+3] - FR[k+3] + FL[k+3] + S[k]);
			UL_star[k+3] = UR_star[k+3];
		}
	}
	if(sens_hllc)
	{
		double s_lambda1 = s_utilde - s_atilde;
		double s_lambda3 = s_utilde + s_atilde;
		double b1 = s_lambda1*(UL_star[0]-UL[0])*d1 + UL[4] - lambda[0]*UL[3];
		double b2 = s_utilde*(UR_star[0] - UL_star[0]);
		double b3 = s_lambda3*(UR[0]-UR_star[0])*d3 - UR[4] + lambda[2]*UR[3];
		double b4 = s_utilde*(UR_star[1]-UL_star[1]);
		double b5 = S[1] - lambda[0]*UL[4] + lambda[2]*UR[4] - FR[4] + FL[4];
		double b6 = S[2] - lambda[0]*UL[5] + lambda[2]*UR[5] - FR[5] + FL[5];
		
		UL_star[3] = ((2*atilde+utilde)*b1 + (atilde+utilde)*b2 + utilde*b3 - b5)/(2*atilde*atilde);
		UL_star[4] = ((utilde*utilde+atilde*utilde)*b1 + (utilde*utilde-atilde*atilde)*b2 + (utilde*utilde-atilde*utilde)*b3 + (atilde-utilde)*b5)/(2*atilde*atilde);
		
		UL_star[5] = ((utilde*utilde*utilde + atilde*utilde*utilde)*(gamma-1)*b1 + (utilde*utilde*utilde*(gamma-1) + (2-gamma)*2*atilde*atilde*utilde)*b2 + (utilde*utilde*utilde-atilde*utilde*utilde)*(gamma-1)*b3 - 2*atilde*atilde*b4 + utilde*utilde*(1-gamma)*b5 + 2*(gamma-1)*atilde*b6)/(4*atilde*atilde*(gamma-1));
		
		UR_star[3] = (b1+b2+b3)/atilde - UL_star[3];
		UR_star[4] = b5/atilde - UL_star[4];
		UR_star[5] = b6/atilde - UL_star[5];
	}
	return transonic_raref;
};

void roe::compute_residual(vector<vector<double> >& R) const
{
	int N = this->get_size();
	R.resize(D);
	for (int k = 0; k < D; ++k)
		R[k].assign(N,0);
	vector<double> UL(D), UR(D), F(D/2), UL_star(D), UR_star(D), U_SL(D), U_SR(D);
	for (int i = 0; i < N+1; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		vector<double> lambda;
		int transonic_raref = this->compute_U_star(UL, UR, UL_star, UR_star, lambda, U_SL, U_SR, i);
		if(CD)
		{
			if(transonic_raref == 0)
			{
				for (int k = 0; k < D; ++k)
				{
					if(i < N)
					{
						R[k][i] += max(lambda[0]-sigma[i], 0.0)*UL[k] + (max(lambda[1] - sigma[i], 0.0)-max(lambda[0] - sigma[i], 0.0))*UL_star[k] + (max(lambda[2] - sigma[i], 0.0)-max(lambda[1] - sigma[i], 0.0))*UR_star[k] - max(lambda[2], sigma[i])*UR[k];
					}
					
					if(i > 0)
					{
						R[k][i-1] += min(lambda[0], sigma[i])*UL[k] + (min(lambda[1] - sigma[i], 0.0) - min(lambda[0] - sigma[i], 0.0))*UL_star[k] + (min(lambda[2]- sigma[i], 0.0)-min(lambda[1] - sigma[i], 0.0))*UR_star[k] - min(lambda[2] - sigma[i], 0.0)*UR[k];
					}
				}
			}
			if(transonic_raref == 1) //left transonic raref
			{
				vector<double> WL, WL_star;
				this->conservative2physical(UL, WL);
				this->conservative2physical(UL_star, WL_star);
				double lambda1L = WL[1] - sqrt(gamma*WL[2]/WL[0]);
				double lambda1L_star = WL_star[1] - sqrt(gamma*WL_star[2]/WL_star[0]);
				for (int k = 0; k < D; ++k)
				{
					if(i < N)
					{
						R[k][i] += lambda1L_star*(U_SL[k] - UL_star[k]) + lambda[1]*(UL_star[k]-UR_star[k]) + lambda[2]*(UR_star[k]-UR[k]);
					}
					if(i > 0)
					{
						R[k][i-1] -= lambda1L*(U_SL[k] - UL[k]);
					}
				}
			}
		}
		else
		{
			if(transonic_raref == 0)
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
			}
			if(transonic_raref == 1) //left transonic raref
			{
				vector<double> WL, WL_star;
				this->conservative2physical(UL, WL);
				this->conservative2physical(UL_star, WL_star);
				double lambda1L = WL[1] - sqrt(gamma*WL[2]/WL[0]);
				double lambda1L_star = WL_star[1] - sqrt(gamma*WL_star[2]/WL_star[0]);
				for (int k = 0; k < D; ++k)
				{
					if(i < N)
					{
						R[k][i] += lambda1L_star*(U_SL[k] - UL_star[k]) + lambda[1]*(UL_star[k]-UR_star[k]) + lambda[2]*(UR_star[k]-UR[k]);
					}
					if(i > 0)
					{
						R[k][i-1] -= lambda1L*(U_SL[k] - UL[k]);
					}
				}
			}
		}
		
		// adding (dx(h))(P-F)
		vector<double> Wi, Ui, Fi;
		this->get_U(Ui,i);
		this->get_W(Wi,i);
		if(i < N)
		{
			double p (Wi[2]), s_p(Wi[5]);
			this->flux(Ui,Fi);
			R[0][i] += -Fi[0]*delta_h[i]/h[i];
			R[1][i] += (p - Fi[1])*delta_h[i]/h[i];
			R[2][i] += -Fi[2]*delta_h[i]/h[i];
			R[3][i] += -Fi[3]*delta_h[i]/h[i] - Fi[0]*(delta_s_h[i]/h[i] - s_h[i]*delta_h[i]/h[i]/h[i]);
			R[4][i] += (s_p - Fi[4])*delta_h[i]/h[i] + (p - Fi[1])*(delta_s_h[i]/h[i] - s_h[i]*delta_h[i]/h[i]/h[i]);
			R[5][i] += -Fi[5]*delta_h[i]/h[i] - Fi[2]*(delta_s_h[i]/h[i] - s_h[i]*delta_h[i]/h[i]/h[i]);
		}
	}
	return;
};
