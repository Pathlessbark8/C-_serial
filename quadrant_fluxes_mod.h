#pragma once
#include<vector>
    	// This module consists of quadrant split fluxes 
    	// with respect to the x-coordinate direction ..
    
    
        void flux_quad_GxI(std::vector<double> G, double nx, double ny, double u1, double u2, double rho, double pr)
    {

    
         
            double  tx, ty, ut, un;
            double  beta;
            double  S1, B1, S2, B2;
            double  A1neg, A2neg;
            double  temp1, temp2, temp3, temp4;
            double  pr_by_rho, u_sqr;
            double  derf;
    
    
            tx = ny;
            ty = -nx;
        
            ut = u1*tx + u2*ty;
            un = u1*nx + u2*ny;


            beta = 0.5*rho/pr;
            S1 = ut*sqrt(beta) ;
            S2 = un*sqrt(beta) ;
            B1 = 0.5*exp(-S1*S1)/sqrt(pi*beta);
            B2 = 0.5*exp(-S2*S2)/sqrt(pi*beta);
            A1neg = 0.5*(1.0 - erf(S1))  ;   
            A2neg = 0.5*(1.0 - erf(S2)) ;    


            pr_by_rho = pr/rho;
            u_sqr = ut*ut + un*un;
    
    	// Expressions for the split fluxes ..	
    
             G[1] = rho*A2neg*(ut*A1neg - B1)  ;
            temp1 = pr_by_rho + ut*ut;
            temp2 = temp1*A1neg - ut*B1;
             G[2] = rho*A2neg*temp2;
    
            temp1 = ut*A1neg - B1;
            temp2 = un*A2neg - B2;
             G[3] = rho*temp1*temp2;
    
            temp1 = (7.0*pr_by_rho) + u_sqr;
            temp2 = 0.5*ut*temp1*A1neg;
     
            temp1 = (6.0*pr_by_rho) + u_sqr;
            temp3 = 0.5*B1*temp1 ;
    
            temp1 = ut*A1neg - B1;
            temp4 = 0.5*rho*un*B2*temp1;
          
             G[4] = rho*A2neg*(temp2 - temp3) - temp4;
        }
    
    
    
        void flux_quad_GxII(std::vector<double>  G, double nx,double ny,double u1,double u2,double rho,double pr)
    {
            double  tx, ty, ut, un;
            double  beta;
            double  S1, B1, S2, B2;
            double  A1pos, A2neg;
            double  temp1, temp2, temp3, temp4;
            double  pr_by_rho, u_sqr;
            double  derf;
    
    
            tx = ny;
            ty = -nx;
    
            ut = u1*tx + u2*ty;
            un = u1*nx + u2*ny;
    
    
            beta = 0.5*rho/pr;
            S1 = ut*sqrt(beta) ;
            S2 = un*sqrt(beta) ;
            B1 = 0.5*exp(-S1*S1)/sqrt(pi*beta);
            B2 = 0.5*exp(-S2*S2)/sqrt(pi*beta);
            A1pos = 0.5*(1.0 + erf(S1))   ;  
            A2neg = 0.5*(1.0 - erf(S2)) ;    
    
            pr_by_rho = pr/rho;
            u_sqr = ut*ut + un*un;
    
    // 	Expressions for the split fluxes ..	
    
             G[1] = rho*A2neg*(ut*A1pos + B1)  ;
    
            temp1 = pr_by_rho + ut*ut;
            temp2 = temp1*A1pos + ut*B1;
             G[2] = rho*A2neg*temp2;
    
            temp1 = ut*A1pos + B1;
            temp2 = un*A2neg - B2;
             G[3] = rho*temp1*temp2;
    
            temp1 = (7.0*pr_by_rho) + u_sqr;
            temp2 = 0.5*ut*temp1*A1pos;
    
            temp1 = (6.0*pr_by_rho) + u_sqr;
            temp3 = 0.5*B1*temp1 ;
    
            temp1 = ut*A1pos + B1;
            temp4 = 0.5*rho*un*B2*temp1;
    
             G[4] = rho*A2neg*(temp2 + temp3) - temp4;
    
    
    }
    
    
    
        void flux_quad_GxIII(std::vector<double>  G, double nx, double ny, double u1, double u2,double rho,double pr)
    
    {
            
            double  tx, ty, ut, un;
            double  beta;
            double  S1, B1, S2, B2;
            double  A1pos, A2pos;
            double  temp1, temp2, temp3, temp4;
            double  pr_by_rho, u_sqr;
            double  derf;
    
    
            tx = ny;
            ty = -nx;
    
            ut = u1*tx + u2*ty;
            un = u1*nx + u2*ny;
    
            beta = 0.5*rho/pr;
            S1 = ut*sqrt(beta) ;
            S2 = un*sqrt(beta) ;
            B1 = 0.5*exp(-S1*S1)/sqrt(pi*beta);
            B2 = 0.5*exp(-S2*S2)/sqrt(pi*beta);
            A1pos = 0.5*(1.0 + erf(S1));     
            A2pos = 0.5*(1.0 + erf(S2)) ;    
    
            pr_by_rho = pr/rho;
            u_sqr = ut*ut + un*un;
    
    	// Expressions for the split fluxes ..	
    
             G[1]= rho*A2pos*(ut*A1pos + B1)  ;
    
            temp1 = pr_by_rho + ut*ut;
            temp2 = temp1*A1pos + ut*B1;
             G[2] = rho*A2pos*temp2;
    
            temp1 = ut*A1pos + B1;
            temp2 = un*A2pos + B2;
             G[3] = rho*temp1*temp2;
    
            temp1 = (7.0*pr_by_rho) + u_sqr;
            temp2 = 0.5*ut*temp1*A1pos;
    
            temp1 = (6.0*pr_by_rho) + u_sqr;
            temp3 = 0.5*B1*temp1 ;
    
            temp1 = ut*A1pos + B1;
            temp4 = 0.5*rho*un*B2*temp1;
    
             G[4] = rho*A2pos*(temp2 + temp3) + temp4;
    
      }
    
    
    
          void flux_quad_GxIV(std::vector<double>  G,double nx,double ny,double u1, double u2,double rho,double pr)
    {
    
            double  tx, ty, ut, un;
            double  beta;
            double  S1, B1, S2, B2;
            double  A1neg, A2pos;
            double  temp1, temp2, temp3, temp4;
            double  pr_by_rho, u_sqr;
            double  derf;
           
          
            tx = ny;
            ty = -nx;
          
            ut = u1*tx + u2*ty;
            un = u1*nx + u2*ny;
          
            beta = 0.5*rho/pr;
            S1 = ut*sqrt(beta) ;
            S2 = un*sqrt(beta) ;
            B1 = 0.5*exp(-S1*S1)/sqrt(pi*beta);
            B2 = 0.5*exp(-S2*S2)/sqrt(pi*beta);
            A1neg = 0.5*(1.0 - erf(S1)) ;    
            A2pos = 0.5*(1.0 + erf(S2)) ;    
          
            pr_by_rho = pr/rho;
            u_sqr = ut*ut + un*un;
          
          	// Expressions for the split fluxes ..	
          
             G[1] = rho*A2pos*(ut*A1neg - B1);  
             
            temp1 = pr_by_rho + ut*ut;
            temp2 = temp1*A1neg - ut*B1;
             G[2] = rho*A2pos*temp2;
          
            temp1 = ut*A1neg - B1;
            temp2 = un*A2pos + B2;
             G[3] = rho*temp1*temp2;
          
            temp1 = (7.0*pr_by_rho) + u_sqr;
            temp2 = 0.5*ut*temp1*A1neg;
          
            temp1 = (6.0*pr_by_rho) + u_sqr;
            temp3 = 0.5*B1*temp1 ;
          
            temp1 = ut*A1neg - B1;
            temp4 = 0.5*rho*un*B2*temp1;
            
             G[4] = rho*A2pos*(temp2 - temp3) + temp4;
           
    }