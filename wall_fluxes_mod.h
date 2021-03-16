#pragma once
#include "data_structure_mod.h"
#include "quadrant_fluxes_mod.h"
#include "split_fluxes_mod.h"
#include "q_variables_mod.h"
#include "limiters_mod.h"
#include <vector>

// CUDA KERNELS

__global__ void wall_dGx_pos_cuda(points &point, double VL_CONST, double gamma_new, int power)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < 1 || i > max_points){
        return;
    }

    if (point.flag_1[i] != 0){
        return;
    }

    double x_i, y_i, x_k, y_k;
    double tx, ty, nx, ny;
    double delx, dely, det, one_by_det;
    double dels, deln;

    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;

    __shared__ double sum_delx_delf[threads_per_block * 4], sum_dely_delf[threads_per_block * 4];
    __shared__ double other_shared[threads_per_block * 8 * 2], qtilde_shared[threads_per_block * 4];

    double dist, weights;
    double dels_weights, deln_weights;


    sum_delx_sqr = 0.0;
    sum_dely_sqr = 0.0;
    sum_delx_dely = 0.0;

    x_i = point.x[i];
    y_i = point.y[i];

    nx = point.nx[i];
    ny = point.ny[i];

    tx = ny;
    ty = -nx;

    for (int j = 0; j < 4; ++j){
        sum_delx_delf[threadIdx.x + blockDim.x * j] = 0;
        sum_dely_delf[threadIdx.x + blockDim.x * j] = 0;
    }

    for (int j = 1; j <= point.xpos_nbhs[i]; j++)
    {
        int k = point.xpos_conn[i][j];

        x_k = point.x[k];
        y_k = point.y[k];

        delx = x_k - x_i;
        dely = y_k - y_i;

        dels = delx * tx + dely * ty;
        deln = delx * nx + dely * ny;

        dist = sqrt(dels * dels + deln * deln);
        weights = pow(dist, power);

        dels_weights = dels * weights;
        deln_weights = deln * weights;

        sum_delx_sqr = sum_delx_sqr + dels * dels_weights;
        sum_dely_sqr = sum_dely_sqr + deln * deln_weights;

        sum_delx_dely = sum_delx_dely + dels * deln_weights;

        for (int r = 0; r < 4; r++){
            other_shared[threadIdx.x + blockDim.x * r] = 0;
        }

        venkat_limiter_cuda(point, other_shared, qtilde_shared, i, VL_CONST, gamma_new, delx, dely);
        flux_quad_GxII_cuda(nx, ny, other_shared, true);

        venkat_limiter_cuda(point, other_shared, qtilde_shared, k, VL_CONST, gamma_new, delx, dely);
        flux_quad_GxII_cuda(nx, ny, other_shared, false);

        for (int r = 0; r < 4; r++)
        {
            sum_delx_delf[threadIdx.x + blockDim.x * r] += other_shared[threadIdx.x + blockDim.x * r] * dels_weights;
            sum_dely_delf[threadIdx.x + blockDim.x * r] += other_shared[threadIdx.x + blockDim.x * r] * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;

    double delta = (2.0 * point.delta[i]);

    for (int j = 0; j < 4; j++)
    {
        point.flux_res[j][i] = delta * (sum_delx_delf[threadIdx.x + blockDim.x * j] * sum_dely_sqr - sum_dely_delf[threadIdx.x + blockDim.x * j] * sum_delx_dely) * one_by_det;
    }
}

__global__ void wall_dGx_neg_cuda(points &point, double VL_CONST, double gamma_new, int power)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < 1 || i > max_points){
        return;
    }

    if (point.flag_1[i] != 0){
        return;
    }

    double x_i, y_i, x_k, y_k;
    double tx, ty, nx, ny;
    double delx, dely, det, one_by_det;
    double dels, deln;

    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;

    __shared__ double sum_delx_delf[threads_per_block * 4], sum_dely_delf[threads_per_block * 4];
    __shared__ double other_shared[threads_per_block * 8 * 2], qtilde_shared[threads_per_block * 4];

    double dist, weights;
    double dels_weights, deln_weights;

    sum_delx_sqr = 0.0;
    sum_dely_sqr = 0.0;
    sum_delx_dely = 0.0;

    x_i = point.x[i];
    y_i = point.y[i];

    nx = point.nx[i];
    ny = point.ny[i];

    tx = ny;
    ty = -nx;

    for (int j = 0; j < 4; ++j){
        sum_delx_delf[threadIdx.x + blockDim.x * j] = 0;
        sum_dely_delf[threadIdx.x + blockDim.x * j] = 0;
    }

    for (int j = 1; j <= point.xneg_nbhs[i]; j++)

    {
        int k = point.xneg_conn[i][j];

        x_k = point.x[k];
        y_k = point.y[k];

        delx = x_k - x_i;
        dely = y_k - y_i;

        dels = delx * tx + dely * ty;
        deln = delx * nx + dely * ny;

        dist = sqrt(dels * dels + deln * deln);
        weights = pow(dist, power);

        dels_weights = dels * weights;
        deln_weights = deln * weights;

        sum_delx_sqr = sum_delx_sqr + dels * dels_weights;
        sum_dely_sqr = sum_dely_sqr + deln * deln_weights;

        sum_delx_dely = sum_delx_dely + dels * deln_weights;

        for (int r = 0; r < 4; r++){
            other_shared[threadIdx.x + blockDim.x * r] = 0;
        }

        venkat_limiter_cuda(point, other_shared, qtilde_shared, i, VL_CONST, gamma_new, delx, dely);
        flux_quad_GxI_cuda(nx, ny, other_shared, true);

        venkat_limiter_cuda(point, other_shared, qtilde_shared, k, VL_CONST, gamma_new, delx, dely);
        flux_quad_GxI_cuda(nx, ny, other_shared, false);

        for (int r = 0; r < 4; r++)
        {
            sum_delx_delf[threadIdx.x + blockDim.x * r] += other_shared[threadIdx.x + blockDim.x * r] * dels_weights;
            sum_dely_delf[threadIdx.x + blockDim.x * r] += other_shared[threadIdx.x + blockDim.x * r] * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;

    double delta = (2.0 * point.delta[i]);

    for (int j = 0; j < 4; j++)
    {
        point.flux_res[j][i] += delta * (sum_delx_delf[threadIdx.x + blockDim.x * j] * sum_dely_sqr - sum_dely_delf[threadIdx.x + blockDim.x * j] * sum_delx_dely) * one_by_det;
    }
}

__global__ void wall_dGy_neg_cuda(points &point, double VL_CONST, double gamma_new, int power)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < 1 || i > max_points){
        return;
    }

    if (point.flag_1[i] != 0){
        return;
    }

    double x_i, y_i, x_k, y_k;
    double tx, ty, nx, ny;
    double delx, dely, det, one_by_det;
    double dels, deln;

    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;

    __shared__ double sum_delx_delf[threads_per_block * 4], sum_dely_delf[threads_per_block * 4];
    __shared__ double other_shared[threads_per_block * 8 * 2], qtilde_shared[threads_per_block * 4];

    double dist, weights;
    double dels_weights, deln_weights;

    sum_delx_sqr = 0.0;
    sum_dely_sqr = 0.0;
    sum_delx_dely = 0.0;

    x_i = point.x[i];
    y_i = point.y[i];

    nx = point.nx[i];
    ny = point.ny[i];

    tx = ny;
    ty = -nx;

    for (int j = 0; j < 4; ++j){
        sum_delx_delf[threadIdx.x + blockDim.x * j] = 0;
        sum_dely_delf[threadIdx.x + blockDim.x * j] = 0;
    }

    for (int j = 1; j <= point.yneg_nbhs[i]; j++)

    {
        int k = point.yneg_conn[i][j];

        x_k = point.x[k];
        y_k = point.y[k];

        delx = x_k - x_i;
        dely = y_k - y_i;

        dels = delx * tx + dely * ty;
        deln = delx * nx + dely * ny;

        dist = sqrt(dels * dels + deln * deln);
        weights = pow(dist, power);

        dels_weights = dels * weights;
        deln_weights = deln * weights;

        sum_delx_sqr = sum_delx_sqr + dels * dels_weights;
        sum_dely_sqr = sum_dely_sqr + deln * deln_weights;

        sum_delx_dely = sum_delx_dely + dels * deln_weights;

        for (int r = 0; r < 4; r++){
            other_shared[threadIdx.x + blockDim.x * r] = 0;
        }

        venkat_limiter_cuda(point, other_shared, qtilde_shared, i, VL_CONST, gamma_new, delx, dely);
        flux_Gyn_cuda(nx, ny, other_shared, true);

        venkat_limiter_cuda(point, other_shared, qtilde_shared, k, VL_CONST, gamma_new, delx, dely);
        flux_Gyn_cuda(nx, ny, other_shared, false);

        for (int r = 0; r < 4; r++)
        {
            sum_delx_delf[threadIdx.x + blockDim.x * r] += other_shared[threadIdx.x + blockDim.x * r] * dels_weights;
            sum_dely_delf[threadIdx.x + blockDim.x * r] += other_shared[threadIdx.x + blockDim.x * r] * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;

    double delta = (2.0 * point.delta[i]);

    for (int j = 0; j < 4; j++)
    {
        point.flux_res[j][i] += delta * (sum_dely_delf[threadIdx.x + blockDim.x * j] * sum_delx_sqr - sum_delx_delf[threadIdx.x + blockDim.x * j] * sum_delx_dely) * one_by_det;
    }
}