#pragma once
#include "data_structure_mod.h"
#include "quadrant_fluxes_mod.h"
#include "split_fluxes_mod.h"
#include "q_variables_mod.h"
#include "limiters_mod.h"
#include <vector>

// This subroutine evaluates the wall flux derivative dGs_pos

void outer_dGx_pos(std::vector<double> G, int i)
{

    int j, k, r;
    double rho, u1, u2, pr;
    double tx, ty, nx, ny;
    double x_i, y_i, x_k, y_k;
    std::vector<double> G_i(5), G_k(5);
    double delx, dely, det, one_by_det;
    double dels, deln;
    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
    std::vector<double> sum_delx_delf(5, 0), sum_dely_delf(5, 0);
    double dist, weights;
    double temp;
    std::vector<double> qtilde_i(5), qtilde_k(5);
    std::vector<double> phi_i(5), phi_k(5);
    double dels_weights, deln_weights;
    std::vector<double> maxi(5), mini(5);

    sum_delx_sqr = 0.0;
    sum_dely_sqr = 0.0;
    sum_delx_dely = 0.0;

    x_i = point.x[i];
    y_i = point.y[i];

    nx = point.nx[i];
    ny = point.ny[i];

    tx = ny;
    ty = -nx;

    for (j = 1; j <= point.xpos_nbhs[i]; j++)
    {
        k = point.xpos_conn[i][j];

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

        for (int j = 1; j <= 4; j++)
        {
            qtilde_i[j] = point.q[j][i] - 0.5 * (delx * point.dq[1][j][i] + dely * point.dq[2][j][i]);
            qtilde_k[j] = point.q[j][k] - 0.5 * (delx * point.dq[1][j][k] + dely * point.dq[2][j][k]);
        }
        venkat_limiter(qtilde_i, phi_i, i);
        venkat_limiter(qtilde_k, phi_k, k);

        for (int j = 1; j <= 4; j++)
        {
            qtilde_i[j] = point.q[j][i] - 0.5 * phi_i[j] * (delx * point.dq[1][j][i] + dely * point.dq[2][j][i]);
            qtilde_k[j] = point.q[j][k] - 0.5 * phi_k[j] * (delx * point.dq[1][j][k] + dely * point.dq[2][j][k]);
        }

        qtilde_to_primitive(qtilde_i, u1, u2, rho, pr);
        flux_quad_GxIII(G_i, nx, ny, u1, u2, rho, pr);

        qtilde_to_primitive(qtilde_k, u1, u2, rho, pr);
        flux_quad_GxIII(G_k, nx, ny, u1, u2, rho, pr);
        for (int j = 1; j <= 4; j++)
        {
            sum_delx_delf[j] = sum_delx_delf[j] + (G_k[j] - G_i[j]) * dels_weights;
            sum_dely_delf[j] = sum_dely_delf[j] + (G_k[j] - G_i[j]) * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;
    for (int j = 1; j <= 4; j++)
    {
        G[j] = (sum_delx_delf[j] * sum_dely_sqr - sum_dely_delf[j] * sum_delx_dely) * one_by_det;
    }
}

// This subroutine evaluates the wall flux derivative dGs_neg

void outer_dGx_neg(std::vector<double> G, int i)
{
    int j, k, r;
    double rho, u1, u2, pr;
    double x_i, y_i, x_k, y_k;
    double tx, ty, nx, ny;
    std::vector<double> G_i(5), G_k(5);
    double delx, dely, det, one_by_det;
    double dels, deln;

    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
    std::vector<double> sum_delx_delf(5, 0), sum_dely_delf(5, 0);
    double dist, weights;
    double temp;
    std::vector<double> qtilde_i(5), qtilde_k(5);
    std::vector<double> phi_i(5), phi_k(5);
    std::vector<double> maxi(5), mini(5);
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

    for (j = 1; j <= point.xneg_nbhs[i]; j++)

    {
        k = point.xneg_conn[i][j];

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

        for (int j = 1; j <= 4; j++)
        {
            qtilde_i[j] = point.q[j][i] - 0.5 * (delx * point.dq[1][j][i] + dely * point.dq[2][j][i]);
            qtilde_k[j] = point.q[j][k] - 0.5 * (delx * point.dq[1][j][k] + dely * point.dq[2][j][k]);
        }
        venkat_limiter(qtilde_i, phi_i, i);
        venkat_limiter(qtilde_k, phi_k, k);

        for (int j = 1; j <= 4; j++)
        {
            qtilde_i[j] = point.q[j][i] - 0.5 * phi_i[j] * (delx * point.dq[1][j][i] + dely * point.dq[2][j][i]);
            qtilde_k[j] = point.q[j][k] - 0.5 * phi_k[j] * (delx * point.dq[1][j][k] + dely * point.dq[2][j][k]);
        }

        qtilde_to_primitive(qtilde_i, u1, u2, rho, pr);
        flux_quad_GxIV(G_i, nx, ny, u1, u2, rho, pr);

        qtilde_to_primitive(qtilde_k, u1, u2, rho, pr);
        flux_quad_GxIV(G_k, nx, ny, u1, u2, rho, pr);
        for (int j = 1; j <= 4; j++)
        {
            sum_delx_delf[j] = sum_delx_delf[j] + (G_k[j] - G_i[j]) * dels_weights;
            sum_dely_delf[j] = sum_dely_delf[j] + (G_k[j] - G_i[j]) * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;

    for (int j = 1; j <= 4; j++)
    {
        G[j] = (sum_delx_delf[j] * sum_dely_sqr - sum_dely_delf[j] * sum_delx_dely) * one_by_det;
    }
}

void outer_dGy_pos(std::vector<double> G, int i)

{

    int j, k, r;
    double rho, u1, u2, pr;
    double x_i, y_i, x_k, y_k;
    double tx, ty, nx, ny;
    std::vector<double> G_i(5), G_k(5);
    double delx, dely, det, one_by_det;
    double dels, deln;

    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
    std::vector<double> sum_delx_delf(5, 0), sum_dely_delf(5, 0);
    double dist, weights;
    double temp;
    std::vector<double> qtilde_i(5), qtilde_k(5);
    std::vector<double> phi_i(5), phi_k(5);
    std::vector<double> maxi(5), mini(5);
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

    for (j = 1; j <= point.ypos_nbhs[i]; j++)

    {
        k = point.ypos_conn[i][j];

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

        for (int j = 1; j <= 4; j++)
        {
            qtilde_i[j] = point.q[j][i] - 0.5 * (delx * point.dq[1][j][i] + dely * point.dq[2][j][i]);
            qtilde_k[j] = point.q[j][k] - 0.5 * (delx * point.dq[1][j][k] + dely * point.dq[2][j][k]);
        }
        venkat_limiter(qtilde_i, phi_i, i);
        venkat_limiter(qtilde_k, phi_k, k);

        for (int j = 1; j <= 4; j++)
        {
            qtilde_i[j] = point.q[j][i] - 0.5 * phi_i[j] * (delx * point.dq[1][j][i] + dely * point.dq[2][j][i]);
            qtilde_k[j] = point.q[j][k] - 0.5 * phi_k[j] * (delx * point.dq[1][j][k] + dely * point.dq[2][j][k]);
        }

        qtilde_to_primitive(qtilde_i, u1, u2, rho, pr);
        flux_Gyp(G_i, nx, ny, u1, u2, rho, pr);

        qtilde_to_primitive(qtilde_k, u1, u2, rho, pr);
        flux_Gyp(G_k, nx, ny, u1, u2, rho, pr);
        for (int j = 1; j <= 4; j++)
        {
            sum_delx_delf[j] = sum_delx_delf[j] + (G_k[j] - G_i[j]) * dels_weights;
            sum_dely_delf[j] = sum_dely_delf[j] + (G_k[j] - G_i[j]) * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;

    for (int j = 1; j <= 4; j++)
    {
        G[j] = (sum_delx_delf[j] * sum_dely_sqr - sum_dely_delf[j] * sum_delx_dely) * one_by_det;
    }
}
