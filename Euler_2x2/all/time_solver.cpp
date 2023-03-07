#include"time_solver.hpp"
#include <algorithm>

double time_solver::solve (state& st, bool only_state)
{
	bool first_time (true);
	int cont(0);
	double t(start_time), dt, dx(this->m.get_Dx());
	vector<double> lambda;
    bool CD_state = st.get_CD_state();
    bool CD_sens = st.get_CD_sens();
	int N(st.get_size()), D(st.get_dimension());
	if (only_state) D = 2;
	vector<vector<double> > U, Uold, U_int, R;
	st.get_U(Uold);
	st.get_U(U);
	st.get_U(U_int);
	
	/***** CD *****/
	vector<double> x_bar(N+1,0), sigma(N+1,0);
	x_bar[0] = this->m.get_xa(); x_bar[N] = this->m.get_xb();
	vector<vector<double> > U_bar(D);
	st.get_U(U_bar);
	/**************/

	// time loop
	while (t < end_time)
	{
		++cont;
		if(cont%1==0)
		{
            st.print_physical(path, ios::out | ios::app);
		}
		/*** dt computation ***/
		st.compute_lambda(lambda);
		dt = dx*cfl/(*max_element(lambda.begin(),lambda.end()));
		if(first_time && (t+dt)>end_time)
		{
			dt = end_time-t;
			first_time = false;
		}
		/**********************/
		if (order == 2)
		{
			if (CD_state && CD_sens)
			{
				/** staggered grid definition **/
				sigma.assign(N+1,0);
				for (int i = 1; i < N; ++i)
				{
					if(U[1][i] - U[1][i-1] < -1e-2) // shock
					{
						if(U[0][i] - U[0][i-1] < -1e-2 ) //1-shock
							sigma[i] = -lambda[i];
						if(U[0][i] - U[0][i-1] > 1e-2) //2-shock
							sigma[i] = lambda[i];
					}
					x_bar[i] = dx*i + sigma[i]*0.5*dt;
				}
				/********************************/
				st.set_sigma(sigma);
				st.compute_residual(R);
				
				/********* compute U_int ********/
				for (int i=0; i<N; ++i)
				{
					double dxi = dx;
					double coeff = 0.0;
					if (fabs(sigma[i]) > 1e-8 || fabs(sigma[i+1]) > 1e-8)
					{
						dxi = (x_bar[i+1]-x_bar[i]);
						coeff = (dx-dxi)/dxi;
					}
					for (int k = 0; k < 4; ++k)
						U_int[k][i] = Uold[k][i] + coeff*Uold[k][i] + 0.5*dt/dxi*R[k][i];
				}
				/********************************/

				/********* compute U_bar ********/
				st.set_U(U_int);
				st.compute_residual(R);

				for (int i=0; i<N; ++i)
					x_bar[i] = dx*i + sigma[i]*dt;
				
				for (int i=0; i<N; ++i)
				{
					double dxi = dx;
					double coeff = 0.0;
					if (fabs(sigma[i]) > 1e-8 || fabs(sigma[i+1]) > 1e-8)
					{
						dxi = (x_bar[i+1]-x_bar[i]);
						coeff = (dx-dxi)/dxi;
					}
					for (int k = 0; k < 4; ++k)
						U_bar[k][i] = Uold[k][i] + coeff*U_int[k][i] + dt/dxi*R[k][i];
				}
				/********************************/

				/*********** sampling ***********/
				double an;
                can(cont, an);

				for (int i=0; i<N; ++i)
				{
                    for (int k=0; k<4; ++k)
                    {
                        if (an < dt/dx*max(0.0, sigma[i]))
                            U[k][i] = U_bar[k][i-1];
                        
                        if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i];
                        
                        if (an > 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i+1];
                    }
				}
				/********************************/
			}
			else
			{
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
				{
					for (int k=0; k<4; ++k)
						U_int[k][i] = Uold[k][i] + 0.5*dt/dx*R[k][i];
					//cout << i << " " << R[0][i] << "\n";
				}
				//cout << "\nsecondo mezzo time step\n";
				st.set_U(U_int);
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
					for (int k=0; k<4; ++k)
						U[k][i] = Uold[k][i] + dt/dx*R[k][i];
			}
		}
		if(order == 1)
		{
			if(CD_state && CD_sens)
			{
				/** staggered grid definition **/
				sigma.assign(N+1,0);
				for (int i = 1; i < N; ++i)
				{
					if(U[1][i] - U[1][i-1] < -1e-3) // shock
					{
						if(U[0][i] - U[0][i-1] < -1e-3) //1-shock
							sigma[i] = -lambda[i];
						
						if(U[0][i] - U[0][i-1] > 1e-3) //2-shock
							sigma[i] = lambda[i];
						
					}
					x_bar[i] = dx*i + sigma[i]*dt;
				}
				/********************************/
				st.set_sigma(sigma);
				st.compute_residual(R);
				/********* compute U_bar ********/
				for (int i=0; i<N; ++i)
				{
					double dxi = (x_bar[i+1]-x_bar[i]);
					for (int k = 0; k < 4; ++k)
						U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
				}
				/********************************/
				
				/*********** sampling ***********/
                double an;
                double an2;
                can(cont, an);
                can(cont, an2);
                an2 = an;
				for (int i=0; i<N; ++i)
				{
                    for (int k=0; k<2; ++k)
                    {
                        if (an < dt/dx*max(0.0, sigma[i]))
                            U[k][i] = U_bar[k][i-1];
                        
                        if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i];
                        
                        if (an > 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i+1];
                    }
                    
                    for (int k=2; k<4; ++k)
                    {
                        if (an2 < dt/dx*max(0.0, sigma[i]))
                            U[k][i] = U_bar[k][i-1];
                        
                        if (an2 > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i];
                        
                        if (an2 > 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i+1];
                    }
				}
				/********************************/
			}
			else
			{
                if(!CD_state && !CD_sens)
                {
                    st.compute_residual(R);
                    for (int i=0; i<N; ++i)
                        for (int k=0; k<4; ++k)
                            U[k][i] = Uold[k][i] + dt/dx*R[k][i];
                }
                if(CD_state && !CD_sens)
                {
                    sigma.assign(N+1,0);
                    for (int i = 1; i < N; ++i)
                    {
                        if(U[1][i] - U[1][i-1] < -1e-3) // shock
                        {
                            if(U[0][i] - U[0][i-1] < -1e-3) //1-shock
                                sigma[i] = -lambda[i];
                            
                            if(U[0][i] - U[0][i-1] > 1e-3) //2-shock
                                sigma[i] = lambda[i];
                            
                        }
                        x_bar[i] = dx*i + sigma[i]*dt;
                    }
                    /********************************/
                    st.set_sigma(sigma);
                    st.compute_residual(R);
                    /********* compute U_bar ********/
                    for (int i=0; i<N; ++i)
                    {
                        double dxi = (x_bar[i+1]-x_bar[i]);
                        for (int k = 0; k < 2; ++k)
                            U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
                    }
                    /********************************/
                    
                    /*********** sampling ***********/
                    double an;
                    can(cont, an);

                    for (int i=0; i<N; ++i)
                    {
                        for (int k=0; k<2; ++k)
                        {
                            if (an < dt/dx*max(0.0, sigma[i]))
                                U[k][i] = U_bar[k][i-1];
                            
                            if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
                                U[k][i] = U_bar[k][i];
                            
                            if (an > 1+dt/dx*min(0.0, sigma[i+1]))
                                U[k][i] = U_bar[k][i+1];
                            
                        }
                        for (int k=2; k<4; ++k)
                            U[k][i] = Uold[k][i] + dt/dx*R[k][i];
                    }
                }
                if(!CD_state && CD_sens)
                {
                    sigma.assign(N+1,0);
                    for (int i = 1; i < N; ++i)
                    {
                        if(U[1][i] - U[1][i-1] < -1e-3) // shock
                        {
                            if(U[0][i] - U[0][i-1] < -1e-3) //1-shock
                                sigma[i] = -lambda[i];
                            
                            if(U[0][i] - U[0][i-1] > 1e-3) //2-shock
                                sigma[i] = lambda[i];
                            
                        }
                        x_bar[i] = dx*i + sigma[i]*dt;
                    }
                    /********************************/
                    st.set_sigma(sigma);
                    st.compute_residual(R);
                    for (int i=0; i<N; ++i)
                        for (int k=0; k<2; ++k)
                            U[k][i] = Uold[k][i] + dt/dx*R[k][i];
                    /********* compute U_bar ********/
                    for (int i=0; i<N; ++i)
                    {
                        double dxi = (x_bar[i+1]-x_bar[i]);
                        for (int k = 2; k < 4; ++k)
                            U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
                    }
                    /********************************/
                    
                    /*********** sampling ***********/
                    double an;
                    can(cont, an);
                    
                    for (int i=0; i<N; ++i)
                    {
                        for (int k=2; k<4; ++k)
                        {
                            if (an < dt/dx*max(0.0, sigma[i]))
                                U[k][i] = U_bar[k][i-1];
                            
                            if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
                                U[k][i] = U_bar[k][i];
                            
                            if (an > 1+dt/dx*min(0.0, sigma[i+1]))
                                U[k][i] = U_bar[k][i+1];
                        }
                    }
                }
			}
		}
		
		Uold = U;
		st.set_U(U);
		t += dt;
	}
	return t;
}
