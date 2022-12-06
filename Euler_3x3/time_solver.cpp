#include"time_solver.hpp"
#include <algorithm>

double time_solver::solve (state& st, bool only_state)
{
	bool first_time (true);
	int cont(0);
	double t(start_time);
	double dt, dx(this->m.get_Dx());
	vector<vector<double> > U, Uold, W, R, S;
	st.get_U(Uold);
	st.get_U(U);
    bool CD_state = st.get_CD_state();
    bool CD_sens = st.get_CD_sens();
	int N(st.get_size()), D(st.get_dimension());
	bool stationary (false);
	vector<double> maxR(D, 0);

	while (t < end_time && !stationary)
	{
	/*	if (cont%100==0)
		{
			cout << t << endl;
			//cerr << t << endl;
			st.print_physical("", ios::out | ios::app);
		}
	*/	maxR.assign(D, 0);

		++cont;
		vector<int> d1, d3, c;
		st.detector_s1(d1, 4*dx);
		st.detector_s3(d3, 4*dx);
		st.detector_c(c, 4*dx);
		
		/***** CD *****/
		vector<double> x_bar(N+1,0), sigma(N+1,0);
		x_bar[0] = this->m.get_xa(); x_bar[N] = this->m.get_xb();
		vector<vector<double> > U_bar(D), U_int(D);
		st.get_U(U_int);
		st.get_U(U_bar);
		/**************/
		
		/*** dt computation ***/
		double mv = st.compute_maxvel();
		dt = dx*cfl/mv;
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
                if (only_state) D = D/2;
                
				/** staggered grid definition **/
				vector<double> sigma;
				sigma.assign(N+1,0);
				vector<double> l1, l2, l3;
				st.compute_lambda1(l1);
				st.compute_lambda2(l2);
				st.compute_lambda3(l3);
				for (int i = 1; i < N; ++i)
				{
					if( c[i] )
						sigma[i] = 0;//l2[i];
					
					if( d1[i] ) //1-shock
						sigma[i] = l1[i];
					
					if( d3[i] ) //3-shock
						sigma[i] = l3[i];
					x_bar[i] = dx*i + sigma[i]*dt*0.5;
				}
				/********************************/
				st.set_sigma(sigma);
				st.compute_residual(R);
				
				/********* compute U_int ********/
				for (int i = 0; i < N; ++i)
				{
					double dxi = dx;
					double coeff = 0.0;
					if (fabs(sigma[i]) > 1e-8 || fabs(sigma[i+1]) > 1e-8)
					{
						dxi = (x_bar[i+1]-x_bar[i]);
						coeff = (dx-dxi)/dxi;
					}
					for (int k = 0; k < D; ++k)
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
					for (int k = 0; k < D; ++k)
						U_bar[k][i] = Uold[k][i] + coeff*U_int[k][i] + dt/dxi*R[k][i];
				}
				/********************************/
				
				/*********** sampling ***********/
				double an;
				can(cont, an);
				
				for (int i=0; i<N; ++i)
				{
					for (int k=0; k<D; ++k)
					{
						if (an < dt/dx*max(0.0, sigma[i]))
							U[k][i] = U_bar[k][i-1];
						
						if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i];
						
						if (an > 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i+1];
						
					}
				}
			} // end if CD_state && CD_sens
			if( !CD_state && !CD_sens)
			{
                if ( only_state ) D = D/2;
                
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
					for (int k=0; k<D; ++k)
						U_int[k][i] = Uold[k][i] + 0.5*dt/dx*R[k][i];
				
				vector<vector<double> > R_old = R;
				st.set_U(U_int);
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
					for (int k=0; k<D; ++k)
						U[k][i] = Uold[k][i] + 0.5*dt/dx*(R[k][i]+R_old[k][i]);
			}
            if(CD_state && !CD_sens)
            {
                // state solution
                /** staggered grid definition **/
                vector<double> sigma;
                sigma.assign(N+1,0);
                vector<double> l1, l2, l3;
                st.compute_lambda1(l1);
                st.compute_lambda2(l2);
                st.compute_lambda3(l3);
                for (int i = 1; i < N; ++i)
                {
                    if( c[i] )
                        sigma[i] = 0;//l2[i];
                    
                    if( d1[i] ) //1-shock
                        sigma[i] = l1[i];
                    
                    if( d3[i] ) //3-shock
                        sigma[i] = l3[i];
                    x_bar[i] = dx*i + sigma[i]*dt*0.5;
                }
                /********************************/
                st.set_sigma(sigma);
                st.compute_residual(R);
                
                /********* compute U_int (Runge Kutta U^n+1/2) ********/
                for (int i = 0; i < N; ++i)
                {
                    double dxi = dx;
                    double coeff = 0.0;
                    if (fabs(sigma[i]) > 1e-8 || fabs(sigma[i+1]) > 1e-8)
                    {
                        dxi = (x_bar[i+1]-x_bar[i]);
                        coeff = (dx-dxi)/dxi;
                    }
                    // state intermediate solution
                    for (int k = 0; k < D/2; ++k)
                        U_int[k][i] = Uold[k][i] + coeff*Uold[k][i] + 0.5*dt/dxi*R[k][i];
                    
                    // sensitivity intermediate solution
                    if(! only_state)
                    {
                        for (int k = D/2; k < D; ++k)
                            U_int[k][i] = Uold[k][i] + 0.5*dt/dx*R[k][i];
                    }
                }
                /********************************/
                
                /********* compute U_bar ********/
                vector<vector<double> > R_old = R;
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
                    for (int k = 0; k < D/2; ++k)
                        U_bar[k][i] = Uold[k][i] + coeff*U_int[k][i] + dt/dxi*R[k][i];
                }
                /********************************/
                
                /*********** sampling ***********/
                double an;
                can(cont, an);
                
                for (int i=0; i<N; ++i)
                {
                    for (int k = 0; k < D/2; ++k)
                    {
                        if (an < dt/dx*max(0.0, sigma[i]))
                            U[k][i] = U_bar[k][i-1];
                        
                        if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i];
                        
                        if (an > 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i+1];
                    }
                }
                if(!only_state)
                {
                    for (int i = 0; i < N; ++i)
                        for (int k = D/2; k < D; ++k)
                            U[k][i] = Uold[k][i] + 0.5*dt/dx*(R[k][i]+R_old[k][i]);
                }
                
            } // end if(CD_state && !CD_sens)
            
            if(!CD_state && CD_sens)
            {
                if(only_state)
                {
                    st.compute_residual(R);
                    for (int i = 0; i < N; ++i)
                        for (int k = 0; k < D/2; ++k)
                            U_int[k][i] = Uold[k][i] + 0.5*dt/dx*R[k][i];
                    vector<vector<double> > R_old = R;
                    st.set_U(U_int);
                    st.compute_residual(R);
                    for (int i = 0; i < N; ++i)
                        for (int k = 0; k < D/2; ++k)
                            U[k][i] = Uold[k][i] + 0.5*dt/dx*(R[k][i]+R_old[k][i]);
                }
                else
                {
                    /** staggered grid definition **/
                    vector<double> sigma;
                    sigma.assign(N+1,0);
                    vector<double> l1, l2, l3;
                    st.compute_lambda1(l1);
                    st.compute_lambda2(l2);
                    st.compute_lambda3(l3);
                    for (int i = 1; i < N; ++i)
                    {
                        if( c[i] )
                            sigma[i] = 0;//l2[i];
                        
                        if( d1[i] ) //1-shock
                            sigma[i] = l1[i];
                        
                        if( d3[i] ) //3-shock
                            sigma[i] = l3[i];
                        x_bar[i] = dx*i + sigma[i]*dt*0.5;
                    }
                    /********************************/
                    st.set_sigma(sigma);
                    st.compute_residual(R);
                    
                    /********* compute U_int ********/
                    for (int i = 0; i < N; ++i)
                    {
                        // state
                        for (int k = 0; k < D/2; ++k)
                            U_int[k][i] = Uold[k][i] + 0.5*dt/dx*R[k][i];
                        
                        double dxi = dx;
                        double coeff = 0.0;
                        if (fabs(sigma[i]) > 1e-8 || fabs(sigma[i+1]) > 1e-8)
                        {
                            dxi = (x_bar[i+1]-x_bar[i]);
                            coeff = (dx-dxi)/dxi;
                        }
                        // sensitivity
                        for (int k = D/2; k < D; ++k)
                            U_int[k][i] = Uold[k][i] + coeff*Uold[k][i] + 0.5*dt/dxi*R[k][i];
                    }
                    /********************************/
                    
                    /********* compute U_bar ********/
                    st.set_U(U_int);
                    vector<vector<double> > R_old = R;
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
                        for (int k = D/2; k < D; ++k)
                            U_bar[k][i] = Uold[k][i] + coeff*U_int[k][i] + dt/dxi*R[k][i];
                    }
                    /********************************/
                    
                    /*********** sampling ***********/
                    double an;
                    can(cont, an);
                    
                    for (int i = 0; i < N; ++i)
                    {
                        // state
                        for (int k = 0; k < D/2; ++k)
                            U[k][i] = Uold[k][i] + 0.5*dt/dx*(R[k][i]+R_old[k][i]);
                        
                        //sensitivity
                        for (int k = D/2 ; k < D; ++k)
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
                    
            } // end if(!CD_state && CD_sens)
            
		} // end if order == 2
		else // if order 1
		{
			if(CD_state && CD_sens)
			{
                if(only_state) D = D/2;

				/** staggered grid definition **/
				sigma.assign(N+1,0);
				vector<double> l1, l2, l3;
				st.compute_lambda1(l1);
				st.compute_lambda2(l2);
				st.compute_lambda3(l3);
				for (int i = 1; i < N; ++i)
				{
					if( c[i] )
						sigma[i] = 0;//l2[i];
					
					if( d1[i] ) //1-shock
						sigma[i] = l1[i];
					
					if( d3[i] ) //3-shock
						sigma[i] = l3[i];
					x_bar[i] = dx*i + sigma[i]*dt;
				}
				/********************************/
				st.set_sigma(sigma);
				st.compute_residual(R);
				/********* compute U_bar ********/
				for (int i=0; i<N; ++i)
				{
					double dxi = (x_bar[i+1]-x_bar[i]);
					for (int k = 0; k < D; ++k)
						U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
				}
				/********************************/

				/*********** sampling ***********/
				double an;
				can(cont, an);
				
				for (int i=0; i<N; ++i)
				{
					for (int k=0; k<D; ++k)
					{
						if (an < dt/dx*max(0.0, sigma[i]))
							U[k][i] = U_bar[k][i-1];
						
						if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i];
						
						if (an > 1+dt/dx*min(0.0, sigma[i+1]))
							U[k][i] = U_bar[k][i+1];
						
					}
				}
			} // end if(CD_state && CD_sens)
			if (!CD_state && !CD_sens)
			{
                if(only_state) D = D/2;
				st.compute_residual(R);
				for (int i=0; i<N; ++i)
					for (int k=0; k<D; ++k)
						U[k][i] = Uold[k][i] + dt/dx*R[k][i];
			} // end if (!CD_state && !CD_sens)
            if (CD_state && !CD_sens )
            {
                /** staggered grid definition **/
                sigma.assign(N+1,0);
                vector<double> l1, l2, l3;
                st.compute_lambda1(l1);
                st.compute_lambda2(l2);
                st.compute_lambda3(l3);
                for (int i = 1; i < N; ++i)
                {
                    if( c[i] )
                        sigma[i] = 0;//l2[i];
                    
                    if( d1[i] ) //1-shock
                        sigma[i] = l1[i];
                    
                    if( d3[i] ) //3-shock
                        sigma[i] = l3[i];
                    x_bar[i] = dx*i + sigma[i]*dt;
                }
                /********************************/
                st.set_sigma(sigma);
                st.compute_residual(R);
                /********* compute U_bar ********/
                for (int i=0; i<N; ++i)
                {
                    double dxi = (x_bar[i+1]-x_bar[i]);
                    for (int k = 0; k < D/2; ++k)
                        U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
                }
                /********************************/

                /*********** sampling ***********/
                double an;
                can(cont, an);
                
                for (int i = 0; i < N; ++i)
                {
                    // state
                    for (int k = 0; k < D/2; ++k)
                    {
                        if (an < dt/dx*max(0.0, sigma[i]))
                            U[k][i] = U_bar[k][i-1];
                        
                        if (an > dt/dx*max(0.0, sigma[i]) && an < 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i];
                        
                        if (an > 1+dt/dx*min(0.0, sigma[i+1]))
                            U[k][i] = U_bar[k][i+1];
                    }
                    // sensitivity
                    if(!only_state)
                    {
                        for (int k = D/2; k < D; ++k)
                            U[k][i] = Uold[k][i] + dt/dx*R[k][i];
                    }
                }
            } // end if (CD_state && !CD_sens )
            
            if (!CD_state && CD_sens )
            {
                if(only_state)
                {
                    st.compute_residual(R);
                    for (int i = 0; i < N; ++i)
                        for (int k = 0; k < D/2; ++k)
                            U[k][i] = Uold[k][i] + dt/dx*R[k][i];
                }
                else
                {
                    /** staggered grid definition **/
                    sigma.assign(N+1,0);
                    vector<double> l1, l2, l3;
                    st.compute_lambda1(l1);
                    st.compute_lambda2(l2);
                    st.compute_lambda3(l3);
                    for (int i = 1; i < N; ++i)
                    {
                        if( c[i] )
                            sigma[i] = 0;//l2[i];
                        
                        if( d1[i] ) //1-shock
                            sigma[i] = l1[i];
                        
                        if( d3[i] ) //3-shock
                            sigma[i] = l3[i];
                        x_bar[i] = dx*i + sigma[i]*dt;
                    }
                    /********************************/
                    st.set_sigma(sigma);
                    st.compute_residual(R);
                    /********* compute U_bar ********/
                    for (int i=0; i<N; ++i)
                    {
                        double dxi = (x_bar[i+1]-x_bar[i]);
                        for (int k = D/2; k < D; ++k)
                            U_bar[k][i] = dx/dxi*Uold[k][i] + dt/dxi*R[k][i];
                    }
                    /********************************/

                    /*********** sampling ***********/
                    double an;
                    can(cont, an);
                    
                    for (int i=0; i<N; ++i)
                    {
                        // state
                        for (int k = 0; k < D/2; ++k)
                            U[k][i] = Uold[k][i] + dt/dx*R[k][i];
                        
                        // sensitivity
                        for (int k = D/2; k < D; ++k)
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
            } // if (! CD_state && CD_sens )
            
		} // end order 1
		Uold = U;
		st.set_U(U);
		t += dt;
		for (int k=0; k<D; ++k)
		{
			for (int i=0; i<N; ++i)
			{
				if(fabs(R[k][i]) > maxR[k])
					maxR[k] = fabs(R[k][i]);
			}
		}
		if(*max_element(maxR.begin(), maxR.end())/dx < 1e-10)
		{
			stationary = true;
		}
	}
	return t;

}
