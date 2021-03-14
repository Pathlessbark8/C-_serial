#pragma once

#include "data_structure_mod.h"
#include "point_normals_mod.h"
#include "generate_connectivity_mod.h"
#include "fpi_solver_mod.h"
#include "initial_conditions_mod.h"
#include <cuda_runtime.h>


void q_lskum()

{
        int i;

        for (i = 1; i <= local_points; i++)
        {
            point.U_old[0][i] = point.prim[0][i];
            point.U_old[1][i] = point.prim[0][i] * point.prim[1][i];
            point.U_old[2][i] = point.prim[0][i] * point.prim[2][i];
            point.U_old[3][i] = 2.5 * point.prim[3][i] + 0.5 * point.prim[0][i] * (point.prim[1][i] * point.prim[1][i] + point.prim[2][i] * point.prim[2][i]);
        }

        //Set U_old to U for first iteration

        points *point_d;
        unsigned long point_size = sizeof(point);
        cudaStream_t stream;
        cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);

        cudaMalloc(&point_d, point_size);
        cudaMemcpy(point_d, &point, point_size, cudaMemcpyHostToDevice);

        cudaDeviceSynchronize();

        fpi_solver_cuda(point_d, stream);

        cudaMemcpy(&point, point_d, point_size, cudaMemcpyDeviceToHost);
        cudaFree(point_d);
}
