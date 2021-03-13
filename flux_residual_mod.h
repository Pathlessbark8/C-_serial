#pragma once

#include "parameter_mod.h"
#include "data_structure_mod.h"
#include "interior_fluxes_mod.h"
#include "wall_fluxes_mod.h"
#include "outer_fluxes_mod.h"
#include <vector>

void cal_flux_residual()

{
    double Gxp[4], Gxn[4], Gyp[4], Gyn[4];

    for (int i = 1; i <= max_points; ++i)
    {
        if (point.flag_1[i] == 0)
        {

            wall_dGx_pos(Gxp, i);

            wall_dGx_neg(Gxn, i);

            wall_dGy_neg(Gyn, i);

            for (int j = 0; j < 4; j++)
            {
                point.flux_res[j][i] = (2.0 * point.delta[i]) * (Gxp[j] + Gxn[j] + Gyn[j]);
            }
        }
        else if (point.flag_1[i] == 1)
        {

            interior_dGx_pos(Gxp, i);

            interior_dGx_neg(Gxn, i);

            interior_dGy_pos(Gyp, i);

            interior_dGy_neg(Gyn, i);

            for (int j = 0; j < 4; j++)
            {
                point.flux_res[j][i] = point.delta[i] * (Gxp[j] + Gxn[j] + Gyp[j] + Gyn[j]);
            }
        }
        else
        {

            outer_dGx_pos(Gxp, i);

            outer_dGx_neg(Gxn, i);

            outer_dGy_pos(Gyp, i);
            for (int j = 0; j < 4; j++)
            {
                point.flux_res[j][i] = point.delta[i] * (Gxp[j] + Gxn[j] + Gyp[j]);
            }
        }
    }
}

__global__ void cal_flux_residual_cuda(points &point, int power, double VL_CONST, double gamma_new)

{
    int bx = blockIdx.x;
    int tx = threadIdx.x;
    int i = bx * blockDim.x + tx;

    if (i < 1 || i > max_points){
        return;
    }

    double Gxp[4], Gxn[4], Gyp[4], Gyn[4];

    if (point.flag_1[i] == 0)
    {

        wall_dGx_pos_cuda(point, Gxp, i, VL_CONST, gamma_new, power);

        wall_dGx_neg_cuda(point, Gxn, i, VL_CONST, gamma_new, power);

        wall_dGy_neg_cuda(point, Gyn, i, VL_CONST, gamma_new, power);

        for (int j = 0; j < 4; j++)
        {
            point.flux_res[j][i] = (2.0 * point.delta[i]) * (Gxp[j] + Gxn[j] + Gyn[j]);
        }
    }
    else if (point.flag_1[i] == 1)
    {

        interior_dGx_pos_cuda(point, Gxp, i, VL_CONST, gamma_new, power);

        interior_dGx_neg_cuda(point, Gxn, i, VL_CONST, gamma_new, power);

        interior_dGy_pos_cuda(point, Gyp, i, VL_CONST, gamma_new, power);

        interior_dGy_neg_cuda(point, Gyn, i, VL_CONST, gamma_new, power);
        
        for (int j = 0; j < 4; j++)
        {
            point.flux_res[j][i] = point.delta[i] * (Gxp[j] + Gxn[j] + Gyp[j] + Gyn[j]);
        }
    }
    else
    {

        outer_dGx_pos_cuda(point, Gxp, i, VL_CONST, gamma_new, power);

        outer_dGx_neg_cuda(point, Gxn, i, VL_CONST, gamma_new, power);

        outer_dGy_pos_cuda(point, Gyp, i, VL_CONST, gamma_new, power);

        for (int j = 0; j < 4; j++)
        {
            point.flux_res[j][i] = point.delta[i] * (Gxp[j] + Gxn[j] + Gyp[j]);
        }
    }
}