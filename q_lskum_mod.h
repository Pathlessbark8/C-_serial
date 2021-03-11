#pragma once

// !	First written on 14.10.2016
// !	updated on Dec 26, 2016
// !	updated on Dec 29, 2016

#include "data_structure_mod.h"
#include "point_normals_mod.h"
#include "generate_connectivity_mod.h"
#include "fpi_solver_mod.h"
#include "initial_conditions_mod.h"

void q_lskum()

{

        int i;

        //Set U_old to U for first iteration
        for (i = 1; i <= local_points; i++)
        {
                point.U_old[1][i] = point.prim[1][i];
                point.U_old[2][i] = point.prim[1][i] * point.prim[2][i];
                point.U_old[3][i] = point.prim[1][i] * point.prim[3][i];
                point.U_old[4][i] = 2.5 * point.prim[4][i] + 0.5 * point.prim[1][i] * (point.prim[2][i] * point.prim[2][i] + point.prim[3][i] * point.prim[3][i]);
        }

        t = 0.0;
        if (restart == 0)
                itr = 0;
        // printf("%d %d\n",itr,max_iters);
        for (it = itr + 1; it <= itr + max_iters; it++)
        {
                // printf("it:%d\n",it);
                fpi_solver(it);
                // printf("%lf\n",residue);
                t = t + dtg;
        }
}
