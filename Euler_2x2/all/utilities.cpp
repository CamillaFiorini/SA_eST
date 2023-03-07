#include "utilities.hpp"
#include<iostream>

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
    /*int bin[31];
       
    an = 0.;
    for (int i = 0; i < 31; i++) bin[i] = 0;
       
    cbin(n,bin);
       
    for (int k = 0; k < 31; k++) { an = an + bin[k]*(pow(2.,0 - k - 1)); }
    */
    an = (double)rand()/(double)RAND_MAX;
}
