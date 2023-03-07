#include"roe_I.hpp"

void roe_I::get_UL_extrapolated (vector<double>& UL, int i) const
{
	UL.resize(4);
	int left = max(i-1, 0);
	for (int k = 0; k < 4; ++k)
	{
		UL[k] = U[k][left];
	}
	return;
};

void roe_I::get_UR_extrapolated (vector<double>& UR, int i) const
{
	UR.resize(4);
	int N = U[0].size();
	int right = min(i, N-1);
	for (int k = 0; k < 4; ++k)
	{
		UR[k] = U[k][right];
	}
	return;
};

void roe_I::compute_residual(vector<vector<double> >& R) const
{
	int N = this->get_size();
	R.resize(4);
	for (int k = 0; k < 4; ++k)
		R[k].assign(N,0);
	
	#pragma omp parallel default(shared)
	{
		int nt = omp_get_num_threads();
		int tid = omp_get_thread_num();
		int ChunkSize = N/nt;
		int CS = ChunkSize;
		if (tid == nt-1)
			ChunkSize = N/nt + N%nt + 1;
		vector<double> to_add(4, 0);
		vector<double> UL(4), UR(4), Ustar(4);
		
		for (int ip = 0; ip < ChunkSize; ++ip)
		{
			int i = ip + tid*CS;
			int left = max(i-1, 0);
			int right = min(i, N-1);
			this->get_U(UL, left);
			this->get_U(UR, right);
			double lambda = this->compute_lambda(UL, UR);
			
			this->compute_U_star(UL, UR, Ustar);
			if (CD_state && CD_sens)
			{
				for (int k=0; k<4; ++k)
				{
					if (i < N)
						R[k][i] += lambda*(Ustar[k]-UR[k]) - sigma[i]*Ustar[k];
					
					if (i > 0)
					{
						if(ip == 0)
							to_add[k] += lambda*(Ustar[k]-UL[k]) + sigma[i]*Ustar[k];
						else
							R[k][i-1] += lambda*(Ustar[k]-UL[k]) + sigma[i]*Ustar[k];
					}
				}
			}
			else
			{
                if(!CD_sens && !CD_state)
                {
                    
                    for (int k=0; k<4; ++k)
                    {
                        if (i < N)
                            R[k][i] += lambda*(Ustar[k]-UR[k]);
                        
                        if (i > 0)
                        {
                            if(ip == 0)
                                to_add[k] += lambda*(Ustar[k]-UL[k]);
                            else
                                R[k][i-1] += lambda*(Ustar[k]-UL[k]);
                        }
                    }
                }
                if(CD_state && !CD_sens)
                {
                    for (int k=0; k<2; ++k)
                    {
                        if (i < N)
                            R[k][i] += lambda*(Ustar[k]-UR[k]) - sigma[i]*Ustar[k];
                        
                        if (i > 0)
                        {
                            if(ip == 0)
                                to_add[k] += lambda*(Ustar[k]-UL[k]) + sigma[i]*Ustar[k];
                            else
                                R[k][i-1] += lambda*(Ustar[k]-UL[k]) + sigma[i]*Ustar[k];
                        }
                    }
                    for (int k=2; k<4; ++k)
                    {
                        if (i < N)
                            R[k][i] += lambda*(Ustar[k]-UR[k]);
                        
                        if (i > 0)
                        {
                            if(ip == 0)
                                to_add[k] += lambda*(Ustar[k]-UL[k]);
                            else
                                R[k][i-1] += lambda*(Ustar[k]-UL[k]);
                        }
                    }
                }
                if(!CD_state && CD_sens)
                {
                    for (int k=0; k<2; ++k)
                    {
                        if (i < N)
                            R[k][i] += lambda*(Ustar[k]-UR[k]);
                        
                        if (i > 0)
                        {
                            if(ip == 0)
                                to_add[k] += lambda*(Ustar[k]-UL[k]);
                            else
                                R[k][i-1] += lambda*(Ustar[k]-UL[k]);
                        }
                    }
                    for (int k=2; k<4; ++k)
                    {
                        if (i < N)
                            R[k][i] += lambda*(Ustar[k]-UR[k]) - sigma[i]*Ustar[k];
                        
                        if (i > 0)
                        {
                            if(ip == 0)
                                to_add[k] += lambda*(Ustar[k]-UL[k]) + sigma[i]*Ustar[k];
                            else
                                R[k][i-1] += lambda*(Ustar[k]-UL[k]) + sigma[i]*Ustar[k];
                        }
                    }
                    
                }
                
			}
		}
		if (tid != 0)
			for (int k = 0; k < 4; ++k)
				R[k][tid*CS-1] += to_add[k];
		
	}
	return;
};
