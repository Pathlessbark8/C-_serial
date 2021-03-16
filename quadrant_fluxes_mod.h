#pragma once

#define pi 3.141592653589793238462643383279502884

__device__ void flux_quad_GxI_cuda(double nx, double ny, double other_shared[], bool add_mode)
{
    double tx, ty, ut, un;
    double beta;
    double S1, B1, S2, B2;
    double A1neg, A2neg;
    double temp1, temp2, temp3, temp4;
    double pr_by_rho, u_sqr;

    double u1 = other_shared[threadIdx.x + blockDim.x * 4];
    double u2 = other_shared[threadIdx.x + blockDim.x * 5];
    double rho = other_shared[threadIdx.x + blockDim.x * 6];
    double pr = other_shared[threadIdx.x + blockDim.x * 7];

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S1 = ut * sqrt(beta);
    S2 = un * sqrt(beta);
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A1neg = 0.5 * (1.0 - erf(S1));
    A2neg = 0.5 * (1.0 - erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2neg * (ut * A1neg - B1)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2neg * (ut * A1neg - B1)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1neg - ut * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2neg * temp2) + other_shared[threadIdx.x + blockDim.x * 1];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2neg * temp2) - other_shared[threadIdx.x + blockDim.x * 1];
    }

    temp1 = ut * A1neg - B1;
    temp2 = un * A2neg - B2;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1neg;

    temp1 = (6.0 * pr_by_rho) + u_sqr;
    temp3 = 0.5 * B1 * temp1;

    temp1 = ut * A1neg - B1;
    temp4 = 0.5 * rho * un * B2 * temp1;

    other_shared[threadIdx.x + blockDim.x * 3] = rho * A2neg * (temp2 - temp3) - temp4;
}

__device__ void flux_quad_GxII_cuda(double nx, double ny, double other_shared[], bool add_mode)
{
    double tx, ty, ut, un;
    double beta;
    double S1, B1, S2, B2;
    double A1pos, A2neg;
    double temp1, temp2, temp3, temp4;
    double pr_by_rho, u_sqr;

    double u1 = other_shared[threadIdx.x + blockDim.x * 4];
    double u2 = other_shared[threadIdx.x + blockDim.x * 5];
    double rho = other_shared[threadIdx.x + blockDim.x * 6];
    double pr = other_shared[threadIdx.x + blockDim.x * 7];

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S1 = ut * sqrt(beta);
    S2 = un * sqrt(beta);
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A1pos = 0.5 * (1.0 + erf(S1));
    A2neg = 0.5 * (1.0 - erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2neg * (ut * A1pos + B1)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2neg * (ut * A1pos + B1)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1pos + ut * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2neg * temp2) + other_shared[threadIdx.x + blockDim.x * 1];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2neg * temp2) - other_shared[threadIdx.x + blockDim.x * 1];
    }

    temp1 = ut * A1pos + B1;
    temp2 = un * A2neg - B2;
    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1pos;

    temp1 = (6.0 * pr_by_rho) + u_sqr;
    temp3 = 0.5 * B1 * temp1;

    temp1 = ut * A1pos + B1;
    temp4 = 0.5 * rho * un * B2 * temp1;

    other_shared[threadIdx.x + blockDim.x * 3] = rho * A2neg * (temp2 + temp3) - temp4;
}

__device__ void flux_quad_GxIII_cuda(double nx, double ny, double other_shared[], bool add_mode)
{
    double tx, ty, ut, un;
    double beta;
    double S1, B1, S2, B2;
    double A1pos, A2pos;
    double temp1, temp2, temp3, temp4;
    double pr_by_rho, u_sqr;

    double u1 = other_shared[threadIdx.x + blockDim.x * 4];
    double u2 = other_shared[threadIdx.x + blockDim.x * 5];
    double rho = other_shared[threadIdx.x + blockDim.x * 6];
    double pr = other_shared[threadIdx.x + blockDim.x * 7];

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S1 = ut * sqrt(beta);
    S2 = un * sqrt(beta);
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A1pos = 0.5 * (1.0 + erf(S1));
    A2pos = 0.5 * (1.0 + erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2pos * (ut * A1pos + B1)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2pos * (ut * A1pos + B1)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1pos + ut * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2pos * temp2) + other_shared[threadIdx.x + blockDim.x * 1];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2pos * temp2) - other_shared[threadIdx.x + blockDim.x * 1];
    }

    temp1 = ut * A1pos + B1;
    temp2 = un * A2pos + B2;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1pos;

    temp1 = (6.0 * pr_by_rho) + u_sqr;
    temp3 = 0.5 * B1 * temp1;

    temp1 = ut * A1pos + B1;
    temp4 = 0.5 * rho * un * B2 * temp1;

    other_shared[threadIdx.x + blockDim.x * 3] = rho * A2pos * (temp2 + temp3) + temp4;
}

__device__ void flux_quad_GxIV_cuda(double nx, double ny, double other_shared[], bool add_mode)
{
    double tx, ty, ut, un;
    double beta;
    double S1, B1, S2, B2;
    double A1neg, A2pos;
    double temp1, temp2, temp3, temp4;
    double pr_by_rho, u_sqr;

    double u1 = other_shared[threadIdx.x + blockDim.x * 4];
    double u2 = other_shared[threadIdx.x + blockDim.x * 5];
    double rho = other_shared[threadIdx.x + blockDim.x * 6];
    double pr = other_shared[threadIdx.x + blockDim.x * 7];

    tx = ny;
    ty = -nx;

    ut = u1 * tx + u2 * ty;
    un = u1 * nx + u2 * ny;

    beta = 0.5 * rho / pr;
    S1 = ut * sqrt(beta);
    S2 = un * sqrt(beta);
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A1neg = 0.5 * (1.0 - erf(S1));
    A2pos = 0.5 * (1.0 + erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2pos * (ut * A1neg - B1)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * A2pos * (ut * A1neg - B1)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1neg - ut * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2pos * temp2) + other_shared[threadIdx.x + blockDim.x * 1];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * A2pos * temp2) - other_shared[threadIdx.x + blockDim.x * 1];
    }

    temp1 = ut * A1neg - B1;
    temp2 = un * A2pos + B2;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1 * temp2) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1neg;

    temp1 = (6.0 * pr_by_rho) + u_sqr;
    temp3 = 0.5 * B1 * temp1;

    temp1 = ut * A1neg - B1;
    temp4 = 0.5 * rho * un * B2 * temp1;

    other_shared[threadIdx.x + blockDim.x * 3] = rho * A2pos * (temp2 - temp3) + temp4;
}