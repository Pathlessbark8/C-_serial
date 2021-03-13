#pragma once

#include "parameter_mod.h"
#include <vector>

void flux_Gxp(double Gxp[], double nx, double ny, double u1, double u2, double rho, double pr)
{
    double tx, ty, ut, un;
    double beta, S1, B1, A1pos;
    double temp1, temp2;
    double pr_by_rho, u_sqr;
    double derf;

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S1 = ut * sqrt(beta);
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    A1pos = 0.5 * (1.0 + erf(S1));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    //  Expressions for the split fluxes ..

    Gxp[0] = rho * (ut * A1pos + B1);

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1pos + ut * B1;
    Gxp[1] = rho * temp2;

    temp1 = ut * un * A1pos + un * B1;
    Gxp[2] = rho * temp1;

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1pos;
    temp1 = (6.0 * pr_by_rho) + u_sqr;
    Gxp[3] = rho * (temp2 + 0.5 * temp1 * B1);
}

void flux_Gxn(double Gxn[], double nx, double ny, double u1, double u2, double rho, double pr)

{

    double tx, ty, ut, un;
    double beta, S1, B1, A1neg;
    double temp1, temp2;
    double pr_by_rho, u_sqr;
    double derf;

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S1 = ut * sqrt(beta);
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    A1neg = 0.5 * (1.0 - erf(S1));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    // 	Expressions for the split fluxes ..

    Gxn[0] = rho * (ut * A1neg - B1);

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1neg - ut * B1;
    Gxn[1] = rho * temp2;

    temp1 = ut * un * A1neg - un * B1;
    Gxn[2] = rho * temp1;

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1neg;
    temp1 = (6.0 * pr_by_rho) + u_sqr;
    Gxn[3] = rho * (temp2 - 0.5 * temp1 * B1);
}

void flux_Gyp(double Gyp[], double nx, double ny, double u1, double u2, double rho, double pr)

{

    double tx, ty, ut, un;
    double beta, S2, B2, A2pos;
    double temp1, temp2;
    double pr_by_rho, u_sqr;
    double derf;

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S2 = un * sqrt(beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A2pos = 0.5 * (1.0 + erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    // Expressions for the split fluxes ..

    Gyp[0] = rho * (un * A2pos + B2);

    temp1 = pr_by_rho + un * un;
    temp2 = temp1 * A2pos + un * B2;
    Gyp[2] = rho * temp2;

    temp1 = ut * un * A2pos + ut * B2;
    Gyp[1] = rho * temp1;

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * un * temp1 * A2pos;
    temp1 = (6.0 * pr_by_rho) + u_sqr;
    Gyp[3] = rho * (temp2 + 0.5 * temp1 * B2);
}

void flux_Gyn(double Gyn[], double nx, double ny, double u1, double u2, double rho, double pr)
{
    double tx, ty, ut, un;
    double beta, S2, B2, A2neg;
    double temp1, temp2;
    double pr_by_rho, u_sqr;
    double derf;

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S2 = un * sqrt(beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A2neg = 0.5 * (1.0 - erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    // Expressions for the split fluxes ..

    Gyn[0] = rho * (un * A2neg - B2);

    temp1 = pr_by_rho + un * un;
    temp2 = temp1 * A2neg - un * B2;
    Gyn[2] = rho * temp2;

    temp1 = ut * un * A2neg - ut * B2;
    Gyn[1] = rho * temp1;

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * un * temp1 * A2neg;
    temp1 = (6.0 * pr_by_rho) + u_sqr;
    Gyn[3] = rho * (temp2 - 0.5 * temp1 * B2);
}
