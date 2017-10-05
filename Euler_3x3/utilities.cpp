//
//  utilities.cpp
//  
//
//  Created by Cami Fiorini on 29/11/16.
//
//

#include "utilities.hpp"


void cmax(const int & n, int & maxa)
{
 for (int k = 0; k < 31; k++)
 {
  if ( (n - pow(2.,k)) >= 0 ) maxa = k;
 }
}

void cbin(const int & n, int bin[])
{
 for (int i = 0; i < 31; i++) bin[i] = 0;
	
 int maxa = 0, nup = n;
 cmax(nup, maxa);
	
 for (int i = 0; i < 31; i++)
 {
  if ( (nup > 0) || (nup < 0) )
  {
   cmax(nup,maxa);
   bin[maxa] = 1;
   nup = nup - int(pow(2.,maxa));
  }
 }
}

void can(const int & n, double & an)
{
 int bin[31];
	
 an = 0.;
 for (int i = 0; i < 31; i++) bin[i] = 0;
	
 cbin(n,bin);
	
 for (int k = 0; k < 31; k++) { an = an + bin[k]*(pow(2.,0 - k - 1)); }
}

/********** Restoring ***************/
/*	t = 0.0114115027477935;
	ifstream if_rho ("rho_int.txt");
	ifstream if_p ("p_int.txt");
	ifstream if_u ("u_int.txt");
	ifstream if_s_rho ("s_rho_int.txt");
	ifstream if_s_p ("s_p_int.txt");
	ifstream if_s_u ("s_u_int.txt");
	double dummy; int k(0);
	while(if_rho >> dummy)
	{
 rho0[k] = dummy;
 ++k;
	}
	k = 0;
	while(if_p >> dummy)
	{
 p0[k] = dummy;
 ++k;
	}
	k = 0;
	while(if_u >> dummy)
	{
 u0[k] = dummy;
 ++k;
	}
	k = 0;
	while(if_s_rho >> dummy)
	{
 s_rho0[k] = dummy;
 ++k;
	}
	k = 0;
	while(if_s_p >> dummy)
	{
 s_p0[k] = dummy;
 ++k;
	}
	k = 0;
	while(if_s_u >> dummy)
	{
 s_u0[k] = dummy;
 ++k;
	}
	cout << "End of restoring\n"; */
/****************************************/