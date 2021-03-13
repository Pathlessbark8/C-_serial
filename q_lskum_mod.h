#pragma once

// !	First written on 14.10.2016
// !	updated on Dec 26, 2016
// !	updated on Dec 29, 2016

#include "data_structure_mod.h"
#include "point_normals_mod.h"
#include "generate_connectivity_mod.h"
#include "fpi_solver_mod.h"
#include "initial_conditions_mod.h"
#include <cuda_runtime.h>


// __global__ void test(points* point_d){

//         (*point_d).x[1] += 1000000;

// }

void q_lskum()

{



        int i;

        points *point_d;
        int point_size = sizeof(point);
        cudaStream_t stream;

        cudaMalloc(&point_d, point_size);
        cudaMemcpy(point_d, &point, point_size, cudaMemcpyHostToDevice);

        // test<<<1,1>>>(point_d);

        cudaMemcpy(&point, point_d, point_size, cudaMemcpyDeviceToHost);
        cudaFree(point_d);


        

        // //Set U_old to U for first iteration
        // for (i = 1; i <= local_points; i++)
        // {
        //         point.U_old[0][i] = point.prim[0][i];
        //         point.U_old[1][i] = point.prim[0][i] * point.prim[1][i];
        //         point.U_old[2][i] = point.prim[0][i] * point.prim[2][i];
        //         point.U_old[3][i] = 2.5 * point.prim[3][i] + 0.5 * point.prim[0][i] * (point.prim[1][i] * point.prim[1][i] + point.prim[2][i] * point.prim[2][i]);
        // }

        // t = 0.0;
        // if (restart == 0)
        // {
        //         itr = 0;
        // }
        // for (it = itr + 1; it <= itr + max_iters; it++)
        // {
        //         fpi_solver(it);
        //         t = t + dtg;
        // }
}
