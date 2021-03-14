#pragma once

#include "data_structure_mod.h"
#include "flux_residual_mod.h"
#include "state_update_mod.h"
#include "q_variables_mod.h"
#include "objective_function_mod.h"
#include "post_processing_mod.h"
#include <iostream>
#include <iomanip>
#include "dump.h"
#include <fstream>

#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

using namespace std;

void fpi_solver(int temp1)

{

        int i, rk;
        for (i = 1; i <= max_points; i++)
        {
                for (int j = 0; j < 4; j++)
                {
                        point.prim_old[j][i] = point.prim[j][i];
                }
        }
        func_delta();
        for (rk = 1; rk <= rks; rk++)
        {
                eval_q_variables();

                eval_q_derivatives();

                for (i = 1; i <= inner_iterations; i++)
                {
                        eval_q_inner_loop();

                        eval_update_innerloop();
                }

                cal_flux_residual();

                state_update(rk);
        }

        std::cout << setprecision(13) << scientific;

        res_new = sqrt(sum_res_sqr) / max_points;
        if (temp1 <= 2 && restart == 0)
        {
                res_old = res_new;
                residue = 0;
        }
        else
        {
                residue = log10(res_new / res_old);
        }

        cout << "Iterations "<< it << " Residue : " << residue << "\n";
}

void fpi_solver_cuda(points* point_d, cudaStream_t stream)

{
    if (restart == 0)
    {
        itr = 0;
    }

    remove("residue");

    dim3 threads(threads_per_block, 1, 1);
    dim3 grid(ceil(max_points + 1 / threads.x), 1, 1);

    double *sum_res_sqr_d = NULL;
    // double *sum_res_sqr_h = new double[max_points];
    
    cudaMalloc(&sum_res_sqr_d, sizeof(double) * (max_points + 1));
    cudaDeviceSynchronize();
    
    for (it = itr + 1; it <= itr + max_iters; it++)
    {
        func_delta_cuda<<<grid, threads>>>(*point_d, CFL);
        for (int rk = 1; rk <= rks; rk++){
            eval_q_variables_cuda<<<grid, threads>>>(*point_d);
            eval_q_derivatives_cuda<<<grid, threads>>>(*point_d, power);
            for (int i = 1; i <= inner_iterations; i++)
            {
                eval_q_inner_loop_cuda<<<grid, threads>>>(*point_d, power);
                eval_update_innerloop_cuda<<<grid, threads>>>(*point_d);
            }
            cal_flux_residual_cuda<<<grid, threads>>>(*point_d, power, VL_CONST, gamma_new);
            state_update_cuda<<<grid, threads>>>(*point_d, rk, euler, mach, theta, sum_res_sqr_d);
            cudaDeviceSynchronize();
            sum_res_sqr = thrust::reduce(thrust::cuda::par.on(stream), sum_res_sqr_d, sum_res_sqr_d + max_points, 0.0);
        }

        res_new = sqrt(sum_res_sqr) / max_points;
        if (it <= 2 && restart == 0)
        {
            res_old = res_new;
            residue = 0;
        }
        else
        {
            residue = log10(res_new / res_old);
        }

        cout << "Iterations "<< it << " Residue : " << setprecision(16) << residue << "\n";
        // ofstream outfile;
        // outfile.open("residue", std::ios_base::app);
        // outfile << it << " " << setprecision(16) << residue << "\n";
    }
    cudaFree(sum_res_sqr_d);
}