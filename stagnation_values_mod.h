#pragma once

#include "data_structure_mod.h"
#include <iostream>
#include <vector>

void stagnation_pressure()

{

    int i, indexMin, indexMax;
    double p0_inf, gammaPower, pMin, pMax, p0, mach_t, angle;
    std::vector<double> prim(5);

    gammaPower = gamma_new / (gamma_new - 1);
    p0_inf = pr_inf * (pow((1 + ((gamma_new - 1) / 2) * mach * mach), gammaPower));

    for (i = 1; i <= max_points; i++)
    {
        for (int j = 1; j <= 4; j++)
        {
            prim[j] = point.prim[j][i];
        }

        angle = sqrt(gamma_new * prim[4] / prim[1]);
        mach_t = sqrt(pow(prim[2], 2) + pow(prim[3], 2)) / angle;
        p0 = prim[4] * (pow((1 + ((gamma_new - 1) / 2) * mach_t * mach_t), gammaPower));
        if (i == 1)
        {
            pMin = p0;
            indexMin = i;
            indexMax = i;
            pMax = p0;
        }
        else if (p0 < pMin)
        {
            pMin = p0;
            indexMin = i;
        }
        else if (p0 > pMax)
        {
            pMax = p0;
            indexMax = i;
        }
    }
    printf("Stagnation values are %f %f %f %f %f %f", pMin, pMax, pMin / p0_inf, pMax / p0_inf, indexMin, indexMax);
}

void objective_function_J()
{

    int i;
    double p0_inf, gammaPower, p0, p0_sum, constant, angle, mach_t;
    std::vector<double> prim(5);
    ;
    double total_p0;

    gammaPower = gamma_new / (gamma_new - 1);
    p0_inf = pr_inf * (pow((1 + ((gamma_new - 1) / 2) * mach * mach), gammaPower));

    constant = 1 / (pow(p0_inf, 2) * plen);
    p0_sum = 0.0;

    for (int i = 1; i <= local_points; i++)
    {
        for (int j = 1; j <= 4; j++)
        {
            prim[j] = point.prim[j][i];
        }
        angle = sqrt(gamma_new * prim[4] / prim[1]);
        mach_t = sqrt(pow(prim[2], 2) + pow(prim[3], 2)) / angle;
        p0 = prim[4] * (pow((1 + ((gamma_new - 1) / 2) * mach_t * mach_t), gammaPower));
        p0_sum = p0_sum + pow((p0_inf - p0), 2);
    }
    total_loss_stagpressure = total_p0 * constant;

    // if(rank == 0)
    // write(*,*) "J: ", total_loss_stagpressure
    // endif
}
