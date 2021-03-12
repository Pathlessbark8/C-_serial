#pragma once

#include "data_structure_mod.h"
#include "flux_residual_mod.h"
#include <vector>
#include <iostream>

std::vector<double> return_column(std::vector<std::vector<double>> a, int k)
{
    std::vector<double> temp;
    for (int i = 1; i <= 4; i++)
    {
        temp.push_back(a[i][k]);
    }
    return temp;
}

void primitive_to_conserved(std::vector<double> prim, double nx, double ny, std::vector<double> &U)
{

    double rho;
    double temp1, temp2;
    // std::cout << prim.size() << " " << U.size() << " " << rho << std::endl;
    rho = prim[0];

    U[0] = rho;
    temp1 = rho * prim[1];
    temp2 = rho * prim[2];
    U[3] = 2.5 * prim[3] + 0.5 * (temp1 * temp1 + temp2 * temp2) / rho;

    U[1] = temp1 * ny - temp2 * nx;
    U[2] = temp1 * nx + temp2 * ny;
}

void conserved_vector_Ubar(std::vector<double> prim, std::vector<double> &Ubar, double nx, double ny)

{

    double u1_inf, u2_inf, u1_inf_rot, u2_inf_rot, e_inf;
    double u1, u2, pr, rho, u1_rot, u2_rot, e;
    double beta, S2, B2_inf, A2n_inf;
    double B2, A2p, temp1, temp2;
    double tx, ty;

    u1_inf = q_inf[1];
    u2_inf = q_inf[2];

    tx = ny;
    ty = -nx;

    u1_inf_rot = u1_inf * tx + u2_inf * ty;
    u2_inf_rot = u1_inf * nx + u2_inf * ny;

    temp1 = (u1_inf_rot * u1_inf_rot + u2_inf_rot * u2_inf_rot);
    e_inf = pr_inf / (rho_inf * (gamma_new - 1.0)) + 0.5 * (temp1);

    beta = (0.5 * rho_inf) / pr_inf;
    S2 = u2_inf_rot * sqrt(beta);
    B2_inf = exp(-S2 * S2) / (2.0 * sqrt(pi * beta));
    A2n_inf = 0.5 * (1.0 - erf(S2));

    rho = prim[0];
    u1 = prim[1];
    u2 = prim[2];
    pr = prim[3];

    u1_rot = u1 * tx + u2 * ty;
    u2_rot = u1 * nx + u2 * ny;

    temp1 = (u1_rot * u1_rot + u2_rot * u2_rot);
    e = pr / (rho * (gamma_new - 1.0)) + 0.5 * (temp1);

    beta = (rho) / (2.0 * pr);
    S2 = u2_rot * sqrt(beta);
    B2 = exp(-S2 * S2) / (2.0 * sqrt(pi * beta));
    A2p = 0.5 * (1.0 + erf(S2));

    Ubar[0] = (rho_inf * A2n_inf) + (rho * A2p);

    Ubar[1] = (rho_inf * u1_inf_rot * A2n_inf) + (rho * u1_rot * A2p);

    temp1 = rho_inf * (u2_inf_rot * A2n_inf - B2_inf);
    temp2 = rho * (u2_rot * A2p + B2);
    Ubar[2] = temp1 + temp2;

    temp1 = (rho_inf * A2n_inf * e_inf - 0.5 * rho_inf * u2_inf_rot * B2_inf);
    temp2 = (rho * A2p * e + 0.5 * rho * u2_rot * B2);

    Ubar[3] = temp1 + temp2;
}

void state_update(int rk)

{

    int i, k, r;
    double delt;
    std::vector<double> U(4), U_old(4);
    double res_sqr;
    double temp;
    double nx, ny;
    double U2_rot, U3_rot;
    const double obt = 1.0 / 3.0;
    const double tbt = 2.0 / 3.0;

    max_res = 0.0;
    sum_res_sqr = 0.0;

    // cout<<"obt "<<obt<<endl;
    // cout<<"tbt "<<tbt<<endl;
    for (i = 1; i <= max_points; i++)
    {
        if (point.flag_1[i] == 0)
        {
            k = i;

            nx = point.nx[k];
            ny = point.ny[k];
            std::vector<double> temp1(4);
            for (int r = 0; r < 4; r++)
            {
                temp1[r] = point.prim[r][k];
            }
            std::vector<double> temp2(4);
            for (int r = 0; r < 4; r++)
            {
                temp2[r] = point.prim_old[r][k];
            }
            primitive_to_conserved(temp1, nx, ny, U);
            primitive_to_conserved(temp2, nx, ny, U_old);

            temp = U[0];

            if (rk != 3)
            {
                for (int j = 0; j < 4; j++)
                {
                    U[j] = U[j] - (0.5 * euler * point.flux_res[j][k]);
                }
            }
            else

            {
                for (int j = 0; j <4; j++)
                {
                    U[j] = tbt * U_old[j] + obt * (U[j] - 0.5 * point.flux_res[j][k]);
                }
            }

            U[2] = 0.;

            U2_rot = U[1];
            U3_rot = U[2];
            U[1] = U2_rot * ny + U3_rot * nx;
            U[2] = U3_rot * ny - U2_rot * nx;

            res_sqr = (U[0] - temp) * (U[0] - temp);

            if (res_sqr > max_res)
            {
                max_res = res_sqr;
                max_res_point = k;
            }

            sum_res_sqr = sum_res_sqr + res_sqr;

            point.prim[0][k] = U[0];
            temp = 1.0 / U[0];
            point.prim[1][k] = U[1] * temp;
            point.prim[2][k] = U[2] * temp;
            point.prim[3][k] = 0.4 * U[3] - (0.2 * temp) * (U[1] * U[1] + U[2] * U[2]);
        }
        else if (point.flag_1[i] == 1)
        {
            k = i;

            nx = point.nx[k];
            ny = point.ny[k];
            std::vector<double> temp1(4);
            for (int r = 0; r < 4; r++)
            {
                temp1[r] = point.prim[r][k];
            }
            std::vector<double> temp2(4);
            for (int r = 0; r < 4; r++)
            {
                temp2[r] = point.prim_old[r][k];
            }
            primitive_to_conserved(temp1, nx, ny, U);
            primitive_to_conserved(temp2, nx, ny, U_old);

            temp = U[0];

            if (rk != 3)
            {
                for (int j =0; j < 4; j++)
                {
                    U[j] = U[j] - (0.5 * euler * point.flux_res[j][k]);
                }
            }
            else

            {
                for (int j = 0; j < 4; j++)
                {
                    U[j] = tbt * U_old[j] + obt * (U[j] - 0.5 * point.flux_res[j][k]);
                }
            }
            
            U2_rot = U[1];
            U3_rot = U[2];
            U[1] = U2_rot * ny + U3_rot * nx;
            U[2] = U3_rot * ny - U2_rot * nx;

            res_sqr = (U[0] - temp) * (U[0] - temp);

            if (res_sqr > max_res)
            {
                max_res = res_sqr;
                max_res_point = k;
            }

            sum_res_sqr = sum_res_sqr + res_sqr;

            point.prim[0][k] = U[0];
            temp = 1.0 / U[0];
            point.prim[1][k] = U[1] * temp;
            point.prim[2][k] = U[2] * temp;

            point.prim[3][k] = 0.4 * U[3] - (0.2 * temp) * (U[1] * U[1] + U[2] * U[2]);
        }
        else
        {
            k = i;

            nx = point.nx[k];
            ny = point.ny[k];
            std::vector<double> temp1(4);
            for (int r = 0; r < 4; r++)
            {
                temp1[r] = point.prim[r][k];
            }
            std::vector<double> temp2(4);
            for (int r = 0; r < 4; r++)
            {
                temp2[r] = point.prim_old[r][k];
            }
            conserved_vector_Ubar(temp1, U, nx, ny);
            conserved_vector_Ubar(temp2, U_old, nx, ny);

            temp = U[0];

            if (rk != 3)
            {
                for (int j = 0; j < 4; j++)
                {
                    U[j] = U[j] - (0.5 * euler * point.flux_res[j][k]);
                }
            }
            else

            {
                for (int j = 0; j < 4; j++)
                {
                    U[j] = tbt * U_old[j] + obt * (U[j] - 0.5 * point.flux_res[j][k]);
                }
            }

            U2_rot = U[1];

            U3_rot = U[2];

            U[1] = U2_rot * ny + U3_rot * nx;

            U[2] = U3_rot * ny - U2_rot * nx;

            point.prim[0][k] = U[0];
            temp = 1.0 / U[0];
            point.prim[1][k] = U[1] * temp;
            point.prim[2][k] = U[2] * temp;
            point.prim[3][k] = 0.4 * U[3] - (0.2 * temp) * (U[1] * U[1] + U[2] * U[2]);
        }
    }
    // for (i = 1; i <= wall_points; i++)
    // {
    //     k = wall_points_index[i];

    //     nx = point.nx[k];
    //     ny = point.ny[k];
    //     primitive_to_conserved(point.prim[k], nx, ny, U);
    //     primitive_to_conserved(point.prim_old[k], nx, ny, U_old);

    //     temp = U[0];

    //     if (rk != 3)
    //     {
    //         for (int j = 1; j <= 4; j++)
    //         {
    //             U[j] = U[j] - (0.5 * euler * point.flux_res[j][k]);
    //         }
    //     }
    //     else

    //     {
    //         for (int j = 1; j <= 4; j++)
    //         {
    //             U[j] = tbt * U_old[j] + obt * (U[j] - 0.5 * point.flux_res[j][k]);
    //         }
    //     }

    //     U[2] = 0.;

    //     U2_rot = U[1];
    //     U3_rot = U[2];
    //     U[1] = U2_rot * ny + U3_rot * nx;
    //     U[2] = U3_rot * ny - U2_rot * nx;

    //     res_sqr = (U[0] - temp) * (U[0] - temp);

    //     if (res_sqr > max_res)
    //     {
    //         max_res = res_sqr;
    //         max_res_point = k;
    //     }

    //     sum_res_sqr = sum_res_sqr + res_sqr;

    //     point.prim[0][k] = U[0];
    //     temp = 1.0 / U[0];
    //     point.prim[1][k] = U[1] * temp;
    //     point.prim[2][k] = U[2] * temp;
    //     point.prim[3][k] = 0.4 * U[3] - (0.2 * temp) * (U[1] * U[1] + U[2] * U[2]);
    // }

    // for (i = 1; i <= outer_points; i++)

    // {
    //     k = outer_points_index[i];

    //     nx = point.nx[k];
    //     ny = point.ny[k];

    //     conserved_vector_Ubar(point.prim[k], U, nx, ny);
    //     conserved_vector_Ubar(point.prim_old[k], U_old, nx, ny);

    //     temp = U[0];

    //     if (rk != 3)
    //     {
    //         for (int j = 1; j <= 4; j++)
    //         {
    //             U[j] = U[j] - (0.5 * euler * point.flux_res[j][k]);
    //         }
    //     }
    //     else
    //     {
    //         for (int j = 1; j <= 4; j++)
    //         {
    //             U[j] = tbt * U_old[j] + obt * (U[j] - 0.5 * point.flux_res[j][k]);
    //         }
    //     }

    //     U2_rot = U[1];

    //     U3_rot = U[2];

    //     U[1] = U2_rot * ny + U3_rot * nx;

    //     U[2] = U3_rot * ny - U2_rot * nx;

    //     point.prim[0][k] = U[0];
    //     temp = 1.0 / U[0];
    //     point.prim[1][k] = U[1] * temp;
    //     point.prim[2][k] = U[2] * temp;
    //     point.prim[3][k] = 0.4 * U[3] - (0.2 * temp) * (U[1] * U[1] + U[2] * U[2]);
    // }

    // for (i = 1; i <= interior_points; i++)
    // {
    //     k = interior_points_index[i];

    //     nx = point.nx[k];
    //     ny = point.ny[k];

    //     primitive_to_conserved(point.prim[k], nx, ny, U);
    //     primitive_to_conserved(point.prim_old[k], nx, ny, U_old);

    //     temp = U[0];

    //     if (rk != 3)
    //     {
    //         for (int j = 1; j <= 4; j++)
    //         {
    //             U[j] = U[j] - (0.5 * euler * point.flux_res[j][k]);
    //         }
    //     }
    //     else
    //     {
    //         for (int j = 1; j <= 4; j++)
    //         {
    //             U[j] = tbt * U_old[j] + obt * (U[j] - 0.5 * point.flux_res[j][k]);
    //         }

    //         U2_rot = U[1];
    //         U3_rot = U[2];
    //         U[1] = U2_rot * ny + U3_rot * nx;
    //         U[2] = U3_rot * ny - U2_rot * nx;

    //         res_sqr = (U[0] - temp) * (U[0] - temp);

    //         if (res_sqr > max_res)
    //         {
    //             max_res = res_sqr;
    //             max_res_point = k;
    //         }

    //         sum_res_sqr = sum_res_sqr + res_sqr;

    //         point.prim[0][k] = U[0];
    //         temp = 1.0 / U[0];
    //         point.prim[0][k] = U[1] * temp;
    //         point.prim[0][k] = U[2] * temp;

    //         point.prim[0][k] = 0.4 * U[3] - (0.2 * temp) * (U[1] * U[1] + U[2] * U[2]);
    //     }
    // }
}

void conserved_to_primitive(std::vector<double> U, std::vector<double> prim)
{

    double temp;

    prim[0] = U[0];

    temp = 1.0 / U[0];

    prim[1] = U[1] * temp;
    prim[2] = U[2] * temp;

    temp = U[3] - (0.5 * temp) * (U[1] * U[1] + U[2] * U[2]);

    prim[3] = 0.4 * temp;
}

// This subroutine computes the delta_t (local time step) at a given point ..

void func_delta()
{

    int i, k, r;
    double delta_t;
    double min_dist, lmin = 1.0;
    double gmin;
    double x_i, x_k, y_k;
    double u1, u2, rho, pr, mod_u;
    double dist;
    double y_i;
    double min_delt;

    for (i = 1; i <= max_points; i++)
    {
        min_delt = 1.0;
        for (r = 1; r <= point.nbhs[i]; r++)
        {
            k = point.conn[i][r];

            rho = point.prim[0][k];
            u1 = point.prim[1][k];
            u2 = point.prim[2][k];
            pr = point.prim[3][k];

            x_i = point.x[i];
            y_i = point.y[i];

            x_k = point.x[k];
            y_k = point.y[k];

            dist = (x_k - x_i) * (x_k - x_i) + (y_k - y_i) * (y_k - y_i);

            dist = sqrt(dist);

            mod_u = sqrt(u1 * u1 + u2 * u2);

            delta_t = dist / (mod_u + 3.0 * sqrt(pr / rho));

            delta_t = CFL * delta_t;

            if (min_delt > delta_t)
            {
                min_delt = delta_t;
            }
        }
        point.delta[i] = min_delt;
    }
}
