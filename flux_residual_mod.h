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

    for (int i = 1; i <= max_points; ++i)
    {
        if (point.flag_1[i] == 0)
        {

            wall_dGx_pos(Gxp, i);

            wall_dGx_neg(Gxn, i);

            wall_dGy_neg(Gyn, i);

            for (int j = 1; j <= 4; j++)
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

            for (int j = 1; j <= 4; j++)
            {
                point.flux_res[j][i] = point.delta[i] * (Gxp[j] + Gxn[j] + Gyp[j] + Gyn[j]);
            }
            if (i == 1)
            {
                cout << "GXP\n";
                for (int r = 1; r <= 4; r++)
                {
                    cout << Gxp[r] << " ";
                }
                cout << endl;
                cout << "GXN\n";
                for (int r = 1; r <= 4; r++)
                {
                    cout << Gxn[r] << " ";
                }
                cout << endl;
                cout << "GYP\n";
                for (int r = 1; r <= 4; r++)
                {
                    cout << Gyp[r] << " ";
                }
                cout << endl;
                cout << "GYN\n";
                for (int r = 1; r <= 4; r++)
                {
                    cout << Gyn[r] << " ";
                }
                cout << endl;
            }
        }
        else
        {

            outer_dGx_pos(Gxp, i);

            outer_dGx_neg(Gxn, i);

            outer_dGy_pos(Gyp, i);
            for (int j = 1; j <= 4; j++)
            {
                point.flux_res[j][i] = point.delta[i] * (Gxp[j] + Gxn[j] + Gyp[j]);
            }
        }
    }
}
