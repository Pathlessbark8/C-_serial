#pragma once
#include "data_structure_mod.h"
#include "quadrant_fluxes_mod.h"
#include "split_fluxes_mod.h"
#include "q_variables_mod.h"
#include "limiters_mod.h"
#include <vector>

// This subroutine evaluates the wall flux derivative dGs_pos

void wall_dGx_pos(std::vector<double> &G, int i)
{

    int j, k, r;
    double rho, u1, u2, pr;
    double tx, ty, nx, ny;
    double x_i, y_i, x_k, y_k;
    std::vector<double> G_i(4), G_k(4);
    double delx, dely, det, one_by_det;
    double dels, deln;
    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
    std::vector<double> sum_delx_delf(5, 0), sum_dely_delf(5, 0);
    double dist, weights;
    double temp;
    std::vector<double> qtilde_i(4), qtilde_k(4);
    std::vector<double> phi_i(4), phi_k(4);
    double dels_weights, deln_weights;
    std::vector<double> maxi(4), mini(4);

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

        for (int r = 0; r < 4; r++)
        {
            qtilde_i[r] = point.q[r][i] - 0.5 * (delx * point.dq[0][r][i] + dely * point.dq[1][r][i]);
            qtilde_k[r] = point.q[r][k] - 0.5 * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k]);
        }
        venkat_limiter(qtilde_i, phi_i, i);
        venkat_limiter(qtilde_k, phi_k, k);

        for (int r = 0; r < 4; r++)
        {
            qtilde_i[r] = point.q[r][i] - 0.5 * phi_i[r] * (delx * point.dq[0][r][i] + dely * point.dq[1][r][i]);
            qtilde_k[r] = point.q[r][k] - 0.5 * phi_k[r] * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k]);
        }

        qtilde_to_primitive(qtilde_i, u1, u2, rho, pr);
        flux_quad_GxII(G_i, nx, ny, u1, u2, rho, pr);

        qtilde_to_primitive(qtilde_k, u1, u2, rho, pr);
        flux_quad_GxII(G_k, nx, ny, u1, u2, rho, pr);
        for (int r = 0; r < 4; r++)
        {
            sum_delx_delf[r] = sum_delx_delf[r] + (G_k[r] - G_i[r]) * dels_weights;
            sum_dely_delf[r] = sum_dely_delf[r] + (G_k[r] - G_i[r]) * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;
    for (int j = 0; j < 4; j++)
    {
        G[j] = (sum_delx_delf[j] * sum_dely_sqr - sum_dely_delf[j] * sum_delx_dely) * one_by_det;
    }
}

// This subroutine evaluates the wall flux derivative dGs_neg

void wall_dGx_neg(std::vector<double> &G, int i)
{
    int j, k, r;
    double rho, u1, u2, pr;
    double x_i, y_i, x_k, y_k;
    double tx, ty, nx, ny;
    std::vector<double> G_i(4), G_k(4);
    double delx, dely, det, one_by_det;
    double dels, deln;

    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
    std::vector<double> sum_delx_delf(5, 0), sum_dely_delf(5, 0);
    double dist, weights;
    double temp;
    std::vector<double> qtilde_i(4), qtilde_k(4);
    std::vector<double> phi_i(4), phi_k(4);
    std::vector<double> maxi(4), mini(4);
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

        for (int r = 0; r < 4; r++)
        {
            qtilde_i[r] = point.q[r][i] - 0.5 * (delx * point.dq[0][r][i] + dely * point.dq[1][r][i]);
            qtilde_k[r] = point.q[r][k] - 0.5 * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k]);
        }
        venkat_limiter(qtilde_i, phi_i, i);
        venkat_limiter(qtilde_k, phi_k, k);

        for (int r = 0; r < 4; r++)
        {
            qtilde_i[r] = point.q[r][i] - 0.5 * phi_i[r] * (delx * point.dq[0][r][i] + dely * point.dq[1][r][i]);
            qtilde_k[r] = point.q[r][k] - 0.5 * phi_k[r] * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k]);
        }

        qtilde_to_primitive(qtilde_i, u1, u2, rho, pr);
        flux_quad_GxI(G_i, nx, ny, u1, u2, rho, pr);

        qtilde_to_primitive(qtilde_k, u1, u2, rho, pr);
        flux_quad_GxI(G_k, nx, ny, u1, u2, rho, pr);
        for (int r = 0; r < 4; r++)
        {
            sum_delx_delf[r] = sum_delx_delf[r] + (G_k[r] - G_i[r]) * dels_weights;
            sum_dely_delf[r] = sum_dely_delf[r] + (G_k[r] - G_i[r]) * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;

    for (int j = 0; j < 4; j++)
    {
        G[j] = (sum_delx_delf[j] * sum_dely_sqr - sum_dely_delf[j] * sum_delx_dely) * one_by_det;
    }
}

void wall_dGy_neg(std::vector<double> &G, int i)
{
    int j, k, r;
    double rho, u1, u2, pr;
    double x_i, y_i, x_k, y_k;
    double tx, ty, nx, ny;
    std::vector<double> G_i(4), G_k(4);
    double delx, dely, det, one_by_det;
    double dels, deln;

    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
    std::vector<double> sum_delx_delf(5, 0), sum_dely_delf(5, 0);
    double dist, weights;
    double temp;
    std::vector<double> qtilde_i(4), qtilde_k(4);
    std::vector<double> phi_i(4), phi_k(4);
    std::vector<double> maxi(4), mini(4);
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

    for (j = 1; j <= point.yneg_nbhs[i]; j++)

    {
        k = point.yneg_conn[i][j];

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

        for (int r = 0; r < 4; r++)
        {
            qtilde_i[r] = point.q[r][i] - 0.5 * (delx * point.dq[0][r][i] + dely * point.dq[1][r][i]);
            qtilde_k[r] = point.q[r][k] - 0.5 * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k]);
        }
        venkat_limiter(qtilde_i, phi_i, i);
        venkat_limiter(qtilde_k, phi_k, k);

        for (int r = 0; r < 4; r++)
        {
            qtilde_i[r] = point.q[r][i] - 0.5 * phi_i[r] * (delx * point.dq[0][r][i] + dely * point.dq[1][r][i]);
            qtilde_k[r] = point.q[r][k] - 0.5 * phi_k[r] * (delx * point.dq[0][r][k] + dely * point.dq[1][r][k]);
        }

        qtilde_to_primitive(qtilde_i, u1, u2, rho, pr);
        flux_Gyn(G_i, nx, ny, u1, u2, rho, pr);

        qtilde_to_primitive(qtilde_k, u1, u2, rho, pr);
        flux_Gyn(G_k, nx, ny, u1, u2, rho, pr);
        for (int r = 0; r < 4; r++)
        {
            sum_delx_delf[r] = sum_delx_delf[r] + (G_k[r] - G_i[r]) * dels_weights;
            sum_dely_delf[r] = sum_dely_delf[r] + (G_k[r] - G_i[r]) * deln_weights;
        }
    }

    det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
    one_by_det = 1.0 / det;

    for (int j = 0; j < 4; j++)
    {
        G[j] = (sum_dely_delf[j]*sum_delx_sqr - sum_delx_delf[j]*sum_delx_dely) * one_by_det;
    }
}