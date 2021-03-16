#pragma once

#include "data_structure_mod.h"
#include <vector>

__device__ void venkat_limiter_cuda(points &point, double other_shared[], double qtilde_shared[], int k, double VL_CONST, double gamma_new, double delx, double dely)
{

    double  del_neg, del_pos;
    double epsi = VL_CONST * point.min_dist[k], num, den, temp;
    epsi = pow(epsi, 3);
    double q[4];
    for (int r = 0; r < 4; r++){
        q[r] = point.q[r][k];
        qtilde_shared[threadIdx.x + blockDim.x * r] = 1;  
    }
 
    for (int r = 0; r < 4; r++)
    {
        del_neg = point.q[r][k] - (0.5 * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k] - q[r]));
        if (abs(del_neg) > 1e-5)
        {
            if (del_neg > 0)
            {
                del_pos = point.qm[0][r][k] - q[r];
            }

            else if (del_neg < 0)
            {
                del_pos = point.qm[1][r][k] - q[r];
            }

            num = (del_pos * del_pos) + (epsi * epsi); // Numerator ..
            num = num * del_neg + 2.0 * del_neg * del_neg * del_pos;

            den = del_pos * del_pos + 2.0 * del_neg * del_neg; // Denominator ..
            den = den + del_neg * del_pos + epsi * epsi;
            den = den * del_neg;

            temp = num / den;

            if (temp < 1)
            {
                qtilde_shared[threadIdx.x + blockDim.x * r] = temp;
            }
        }
    }
    for (int r = 0; r < 4; ++r){
        qtilde_shared[threadIdx.x + blockDim.x * r] = point.q[r][k] - 0.5 * (qtilde_shared[threadIdx.x + blockDim.x * r] * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k]));
    }

    double beta = - qtilde_shared[threadIdx.x + blockDim.x * 3] * 0.5;

    temp = 0.5 / beta;

    double u1 = qtilde_shared[threadIdx.x + blockDim.x * 1] * temp;
    double u2 = qtilde_shared[threadIdx.x + blockDim.x * 2] * temp;

    double temp1 = qtilde_shared[threadIdx.x + blockDim.x * 0] + beta * ( u1 * u1 + u2 * u2 );
    double temp2 = temp1 - (log(beta)/(gamma_new-1));
    double rho = exp(temp2);
    double pr = rho * temp;

    other_shared[threadIdx.x + blockDim.x * 4] = u1;
    other_shared[threadIdx.x + blockDim.x * 5] = u2;
    other_shared[threadIdx.x + blockDim.x * 6] = rho;
    other_shared[threadIdx.x + blockDim.x * 7] = pr;
}
