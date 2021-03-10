#pragma once

#include "data_structure_mod.h"

void compute_entropy()
{

        int k;
        double temp1, temp2;
        double total_entropy;

        total_entropy = 0.;

        temp2 = log(pr_inf);

        for (k = 1; k <= max_points; k++)
        {
                temp1 = pow((point.prim[1][k]), gamma_new);
                temp1 = point.prim[4][k] / temp1;
                temp1 = log(temp1);
                point.entropy[k] = abs(temp1 - temp2);

                total_entropy = total_entropy + abs(temp1 - temp2);
        }
}