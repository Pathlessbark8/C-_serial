#pragma once

#include "parameter_mod.h"
#include "data_structure_mod.h"
#include "interior_fluxes_mod.h"
#include "wall_fluxes_mod.h"
#include "outer_fluxes_mod.h"
#include <vector>

void cal_flux_residual()

{
    int i, k;
    std::vector<double> Gxp(5), Gxn(5), Gyp(5), Gyn(5);

    for (int i = 1; i <= wall_points; i++)
    {

        k = wall_points_index[i];

        wall_dGx_pos(Gxp, k);

        wall_dGx_neg(Gxn, k);

        wall_dGy_neg(Gyn, k);

        for (int j = 1; j <= 4; j++)
        {
            point.flux_res[j][k] = Gxp[j] + Gxn[j] + Gyn[j];
        }

        for (int j = 1; j <= 4; j++)
        {
            point.flux_res[j][k] = 2.0 * point.delta[k] * point.flux_res[j][k];
        }
    }

    for (i = 1; i <= outer_points; i++)
    {

        k = outer_points_index[i];

        outer_dGx_pos(Gxp, k);

        outer_dGx_neg(Gxn, k);

        outer_dGy_pos(Gyp, k);
        for (int j = 1; j <= 4; j++)
        {
            point.flux_res[j][k] = point.delta[k] * (Gxp[j] + Gxn[j] + Gyp[j]);
        }
    }

    for (i = 1; i <= interior_points; i++)
    {
        k = interior_points_index[i];

        interior_dGx_pos(Gxp, k);

        interior_dGx_neg(Gxn, k);

        interior_dGy_pos(Gyp, k);

        interior_dGy_neg(Gyn, k);

        for (int j = 1; j <= 4; j++)
        {
            point.flux_res[j][k] = point.delta[k] * (Gxp[j] + Gxn[j] + Gyp[j] + Gyn[j]);
        }
    }
}
