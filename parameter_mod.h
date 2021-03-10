#pragma once
#include<math.h>

const double mach = 0.63;
const double aoa = 2.0;
double theta;
const double rho_inf = 1.0;
const double pr_inf = 1.0/1.4;
const double gamma_new = 1.4;
const double pi=4*atan(1.0);
double q_init[4], q_inf[4];

  void setup_case_parameters()
  {
        
        // calculate theta
        theta = aoa*pi/180;
        
        // Setup initial conditions
        q_init[1] = rho_inf;
        q_init[2] = Mach*cos(theta);
        q_init[3] = Mach*sin(theta);
        q_init[4] = pr_inf;
        
        // Setup free stream conditions
        q_inf[1] = rho_inf;
        q_inf[2] = Mach*cos(theta);
        q_inf[3] = Mach*sin(theta);
        q_inf[4] = pr_inf;
  }