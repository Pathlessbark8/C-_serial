#pragma once

#define pi 3.141592653589793238462643383279502884

__device__ void flux_Gxp_cuda(double nx, double ny, double other_shared[], bool add_mode)
{
    double tx, ty, ut, un;
    double beta, S1, B1, A1pos;
    double temp1, temp2;
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
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    A1pos = 0.5 * (1.0 + erf(S1));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (ut * A1pos + B1)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (ut * A1pos + B1)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1pos + ut * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp2) + other_shared[threadIdx.x + blockDim.x * 1];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp2) - other_shared[threadIdx.x + blockDim.x * 1];
    }

    temp1 = ut * un * A1pos + un * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1pos;
    temp1 = (6.0 * pr_by_rho) + u_sqr;
    
    other_shared[threadIdx.x + blockDim.x * 3] = rho * (temp2 + 0.5 * temp1 * B1);
}

__device__ void flux_Gxn_cuda(double nx, double ny, double other_shared[], bool add_mode)

{

    double tx, ty, ut, un;
    double beta, S1, B1, A1neg;
    double temp1, temp2;
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
    B1 = 0.5 * exp(-S1 * S1) / sqrt(pi * beta);
    A1neg = 0.5 * (1.0 - erf(S1));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (ut * A1neg - B1)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (ut * A1neg - B1)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + ut * ut;
    temp2 = temp1 * A1neg - ut * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp2) + other_shared[threadIdx.x + blockDim.x * 1];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp2) - other_shared[threadIdx.x + blockDim.x * 1];
    }

    temp1 = ut * un * A1neg - un * B1;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp1) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * ut * temp1 * A1neg;
    temp1 = (6.0 * pr_by_rho) + u_sqr;

    other_shared[threadIdx.x + blockDim.x * 3] = rho * (temp2 - 0.5 * temp1 * B1);
}

__device__ void flux_Gyp_cuda(double nx, double ny, double other_shared[], bool add_mode)

{

    double tx, ty, ut, un;
    double beta, S2, B2, A2pos;
    double temp1, temp2;
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
    S2 = un * sqrt(beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A2pos = 0.5 * (1.0 + erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (un * A2pos + B2)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (un * A2pos + B2)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + un * un;
    temp2 = temp1 * A2pos + un * B2;

    temp1 = ut * un * A2pos + ut * B2;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp1) + other_shared[threadIdx.x + blockDim.x * 1];
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp2) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp1) - other_shared[threadIdx.x + blockDim.x * 1];
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp2) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * un * temp1 * A2pos;
    temp1 = (6.0 * pr_by_rho) + u_sqr;

    other_shared[threadIdx.x + blockDim.x * 3] = rho * (temp2 + 0.5 * temp1 * B2);
}

__device__ void flux_Gyn_cuda(double nx, double ny, double other_shared[], bool add_mode)
{
    double tx, ty, ut, un;
    double beta, S2, B2, A2neg;
    double temp1, temp2;
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
    S2 = un * sqrt(beta);
    B2 = 0.5 * exp(-S2 * S2) / sqrt(pi * beta);
    A2neg = 0.5 * (1.0 - erf(S2));

    pr_by_rho = pr / rho;
    u_sqr = ut * ut + un * un;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (un * A2neg - B2)) + other_shared[threadIdx.x + blockDim.x * 0];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 0] = (rho * (un * A2neg - B2)) - other_shared[threadIdx.x + blockDim.x * 0];
    }

    temp1 = pr_by_rho + un * un;
    temp2 = temp1 * A2neg - un * B2;

    temp1 = ut * un * A2neg - ut * B2;

    if (add_mode){
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp1) + other_shared[threadIdx.x + blockDim.x * 1];
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp2) + other_shared[threadIdx.x + blockDim.x * 2];
    }
    else{
        other_shared[threadIdx.x + blockDim.x * 1] = (rho * temp1) - other_shared[threadIdx.x + blockDim.x * 1];
        other_shared[threadIdx.x + blockDim.x * 2] = (rho * temp2) - other_shared[threadIdx.x + blockDim.x * 2];
    }

    temp1 = (7.0 * pr_by_rho) + u_sqr;
    temp2 = 0.5 * un * temp1 * A2neg;
    temp1 = (6.0 * pr_by_rho) + u_sqr;

    other_shared[threadIdx.x + blockDim.x * 3] = rho * (temp2 - 0.5 * temp1 * B2);
}