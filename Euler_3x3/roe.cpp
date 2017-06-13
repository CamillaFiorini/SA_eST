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
	double utilde = (UL[1]/sqrt(UL[0]) + UR[1]/sqrt(UR[0]))/( sqrt(UL[0]) + sqrt(UR[0]) );
	double HL = UL[2]/UL[0] + ((gamma-1)*(UL[2] - 0.5*UL[1]*UL[1]/UL[0]))/UL[0];
	double HR = UR[2]/UR[0] + ((gamma-1)*(UR[2] - 0.5*UR[1]*UR[1]/UR[0]))/UR[0];
	double Htilde = (sqrt(UL[0])*HL + sqrt(UR[0])*HR)/(sqrt(UL[0]) + sqrt(UR[0]));
	double atilde = sqrt((gamma-1)*(Htilde - 0.5*utilde*utilde));
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
	return (UL[1]/sqrt(UL[0]) + UR[1]/sqrt(UR[0]))/( sqrt(UL[0]) + sqrt(UR[0]) );
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
	double utilde = (UL[1]/sqrt(UL[0]) + UR[1]/sqrt(UR[0]))/( sqrt(UL[0]) + sqrt(UR[0]) );
	double HL = UL[2]/UL[0] + ((gamma-1)*(UL[2] - 0.5*UL[1]*UL[1]/UL[0]))/UL[0];
	double HR = UR[2]/UR[0] + ((gamma-1)*(UR[2] - 0.5*UR[1]*UR[1]/UR[0]))/UR[0];
	double Htilde = (sqrt(UL[0])*HL + sqrt(UR[0])*HR)/(sqrt(UL[0]) + sqrt(UR[0]));
	double atilde = sqrt((gamma-1)*(Htilde - 0.5*utilde*utilde));
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
	return 0.5*UR[3]*uR/sqrt(UR[0])/(sqrt(UR[0])+sqrt(UL[0])) + 0.5*UL[3]*uL/sqrt(UL[0])/(sqrt(UR[0])+sqrt(UL[0])) + (sqrt(UR[0])*s_uR + sqrt(UL[0])*s_uL)/(sqrt(UR[0])+sqrt(UL[0])) + u_tilde/(sqrt(UR[0])+sqrt(UL[0]))*(0.5*UR[3]/sqrt(UR[0]) + 0.5*UL[3]/sqrt(UL[0]));
};


double roe::compute_H(const vector<double>& V) const
{
	return V[2]/V[0] + ((gamma-1)*(V[2] - 0.5*V[1]*V[1]/V[0]))/V[0];
};

double roe::compute_s_H(const vector<double>& V) const
{
	double H = this->compute_H(V);
	double s_u = V[4]/V[0] - V[3]*V[1]/(V[0]*V[0]);
	double s_p = (gamma-1)*(V[5] - 0.5*V[3]*V[1]*V[1]/V[0]/V[0] - V[1]*s_u);
	return (V[5] + s_p)/V[0] + V[3]/V[0]*H;
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
	return 0.5*UR[3]*HR/sqrt(UR[0])/(sqrt(UR[0])+sqrt(UL[0])) + 0.5*UL[3]*HL/sqrt(UL[0])/(sqrt(UR[0])+sqrt(UL[0])) + (sqrt(UR[0])*s_HR + sqrt(UL[0])*s_HL)/(sqrt(UR[0])+sqrt(UL[0])) + H_tilde/(sqrt(UR[0])+sqrt(UL[0]))*(0.5*UR[3]/sqrt(UR[0]) + 0.5*UL[3]/sqrt(UL[0]));
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
	
/*	alpha_tilde[0] = ((Htilde*utilde+0.5*utilde*utilde*atilde - utilde*utilde*utilde)*(UR[0] - UL[0]) + (0.5*utilde*utilde-Htilde-utilde*atilde)*(UR[1] - UL[1]) + (atilde)*(UR[2] - UL[2]))/(atilde*(2*Htilde-utilde*utilde));
	alpha_tilde[1] = ((2*(Htilde-utilde*utilde))*(UR[0] - UL[0]) + (2*utilde)*(UR[1] - UL[1]) + (-2)*(UR[2] - UL[2]))/(2*Htilde-utilde*utilde);
	alpha_tilde[2] = ((utilde*(0.5*utilde*utilde + 0.5*atilde*utilde - Htilde))*(UR[0] - UL[0]) + (-(0.5*utilde*utilde - Htilde + utilde*atilde))*(UR[1] - UL[1]) + (atilde)*(UR[2] - UL[2]))/(atilde*(2*Htilde-utilde*utilde));
*/
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

int roe::detector_s2(const vector<double>& UL, const vector<double>& UR, double threshold) const
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

void roe::compute_flux(const vector<double>& UL, const vector<double>& UR, vector<double>& F, vector<double>& s_Ustar) const
{
	vector<double> FL(D), FR(D);
	F.resize(D/2);
	s_Ustar.resize(D/2);
	this->flux(UL, FL);
	this->flux(UR, FR);
	double lambda1 = this->compute_lambda1(UL, UR);
	double lambda2 = this->compute_lambda2(UL, UR);
	double lambda3 = this->compute_lambda3(UL, UR);
	double utilde = this->compute_utilde(UL, UR);
	double atilde = this->compute_atilde(UL, UR);
	double Htilde = this->compute_Htilde(UL, UR);
	vector<double> r1(D), r2(D), r3(D), alpha_tilde(D), S(D/2), UR_star(D/2), UL_star(D/2);
	int d1(this->detector_s1(UL, UR, 1e-3)), d3(this->detector_s2(UL, UR, 1e-3));
	double s_atilde = this->compute_s_atilde(UL, UR);
	double s_utilde = this->compute_s_utilde(UL, UR);
	double s_Htilde = this->compute_s_Htilde(UL, UR);
	
	r1[0] = 1;
	r1[1] = utilde - atilde;
	r1[2] = Htilde - utilde*atilde;
	r1[3] = 0;
	r1[4] = s_utilde - s_atilde;
	r1[5] = s_Htilde - s_utilde*atilde - utilde*s_atilde;
	
	r2[0] = 1;
	r2[1] = utilde;
	r2[2] = 0.5*utilde*utilde;
	r2[3] = 0;
	r2[4] = s_utilde;
	r2[5] = utilde*s_utilde;
	
	r3[0] = 1;
	r3[1] = utilde + atilde;
	r3[2] = Htilde + utilde*atilde;
	r3[3] = 0;
	r3[4] = s_utilde + s_atilde;
	r3[5] = s_Htilde + s_utilde*atilde + utilde*s_atilde;
	
	this->compute_alpha_tilde(UL, UR, alpha_tilde);
	double dummy1, dummy2;
	for (int k = 0; k < D/2; ++k)
	{
		UL_star[k] = UL[k] + alpha_tilde[0]*r1[k];
		UR_star[k] = UR[k] - alpha_tilde[2]*r3[k];
		S[k] = (s_utilde - s_atilde)*(UL[k] - UL_star[k])*d1 + s_utilde*(UL_star[k] - UR_star[k]) + (s_utilde + s_atilde)*(UR_star[k] - UR[k])*d3;
		F[k] = 0.5*(FL[k] + FR[k]) - 0.5*(alpha_tilde[0]*fabs(lambda1)*r1[k] + alpha_tilde[1]*fabs(lambda2)*r2[k] + alpha_tilde[2]*fabs(lambda3)*r3[k]);
		dummy1 = UL[k+3] + alpha_tilde[0]*r1[k+3] + alpha_tilde[3]*r1[k];
		dummy2 = UR[k+3] - alpha_tilde[2]*r3[k+3] - alpha_tilde[5]*r3[k];
		s_Ustar[k] = 0.5*(dummy1+dummy2)/*1.0/(lambda3-lambda1)*(lambda3*UR[k+3] - lambda1*UL[k+3] - FR[k+3] + FL[k+3] + S[k])*/;
	}
	return;
};

void roe::compute_residual(vector<vector<double> >& R) const
{
	int N = this->get_size();
	R.resize(D);
	for (int k = 0; k < D; ++k)
		R[k].assign(N,0);
	vector<double> UL(D), UR(D), F(D/2), s_Ustar(D/2);
	for (int i = 0; i < N+1; ++i)
	{
		this->get_UL_extrapolated(UL, i);
		this->get_UR_extrapolated(UR, i);
		this->compute_flux(UL, UR, F, s_Ustar);
		double lambda1 = this->compute_lambda1(UL, UR);
		double lambda3 = this->compute_lambda3(UL, UR);
		for (int k = 0; k < D/2; ++k)
		{
			if(i < N)
			{
				R[k][i] += F[k];
				/*if (lambda1 < 0 && lambda3 > 0)
					R[k+3][i] += ;
				if(lambda1 > 0)
					R[k+3][i] += ;*/
				R[k+3][i] += max(lambda1, 0.0)*(UL[k+3] - s_Ustar[k]) + max(lambda3,0.0)*(s_Ustar[k]-UR[k+3]);
			}
			
			
			if(i > 0)
			{
				R[k][i-1] -= F[k];
				/*if(lambda3 < 0)
					R[k+3][i-1] -= ;
				if (lambda3 > 0 && lambda1 < 0)
					R[k+3][i-1] -= ;*/
					
				R[k+3][i-1] += min(lambda1, 0.0)*(UL[k+3] - s_Ustar[k]) + min(lambda3,0.0)*(s_Ustar[k]-UR[k+3]);
			}
		}
	}
	return;
};