#pragma once
#include "data_structure_mod.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

void eval_q_variables()
{

    int k;
    double rho, u1, u2, pr, beta;
    double two_times_beta;

    for (int k = 1; k <= max_points; k++)
    {
        rho = point.prim[1][k];
        u1 = point.prim[2][k];
        u2 = point.prim[3][k];
        pr = point.prim[4][k];

        beta = 0.5 * rho / pr;

        point.q[1][k] = log(rho) + (log(beta) * 2.5) - beta * (u1 * u1 + u2 * u2);

        two_times_beta = 2.0 * beta;

        point.q[2][k] = two_times_beta * u1;
        point.q[3][k] = two_times_beta * u2;

        point.q[4][k] = -two_times_beta;
    }
}

void eval_q_derivatives()
{

    int i, k, r, nbh;

    double x_i, y_i, x_k, y_k;

    double delx, dely, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;

    double det, delq, temp;
    double one_by_det;
    fstream ff;
    ff.open("nbhs", ios::out);
    ff << setprecision(14) << std::scientific;
    cout << setprecision(14) << std::scientific;
    for (int i = 1; i <= max_points; i++)
    {
        x_i = point.x[i];
        y_i = point.y[i];
        sum_delx_sqr = 0.;
        sum_dely_sqr = 0.;
        sum_delx_dely = 0.;

        for (int k = 1; k <= 4; k++)
        {
            point.qm[1][k][i] = point.q[k][i];
            point.qm[2][k][i] = point.q[k][i];
        }
        std::vector<double> sum_delx_delq(5, 0.0), sum_dely_delq(5, 0.0);
        for (int k = 1; k <= point.nbhs[i]; k++)
        {
            nbh = point.conn[i][k];
            for (int r = 1; r <= 4; r++)
            {
                if (point.q[r][nbh] > point.qm[1][r][i])
                {
                    point.qm[1][r][i] = point.q[r][nbh];
                }
                if (point.q[r][nbh] < point.qm[2][r][i])
                {
                    point.qm[2][r][i] = point.q[r][nbh];
                }
            }
            x_k = point.x[nbh];
            y_k = point.y[nbh];

            delx = x_k - x_i;
            dely = y_k - y_i;

            dist = sqrt(delx * delx + dely * dely);
            weights = pow(dist, power);

            sum_delx_sqr = sum_delx_sqr + delx * delx * weights;
            sum_dely_sqr = sum_dely_sqr + dely * dely * weights;

            sum_delx_dely = sum_delx_dely + delx * dely * weights;

            for (int r = 1; r <= 4; r++)
            {
                sum_delx_delq[r] = sum_delx_delq[r] + weights * delx * (point.q[r][nbh] - point.q[r][i]);
                sum_dely_delq[r] = sum_dely_delq[r] + weights * dely * (point.q[r][nbh] - point.q[r][i]);
            }
            // if (i == 2)
            // {
            //     for (int r = 1; r <= 4; r++)
            //     {
            //         cout << sum_delx_delq[r] << " ";
            //     }
            //     cout << endl;
            // }
        }

        det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
        one_by_det = 1.0 / det;

        for (int k = 1; k <= 4; k++)
        {
            point.dq[1][k][i] = (sum_delx_delq[k] * sum_dely_sqr - sum_dely_delq[k] * sum_delx_dely) * one_by_det;
            point.dq[2][k][i] = (sum_dely_delq[k] * sum_delx_sqr - sum_delx_delq[k] * sum_delx_dely) * one_by_det;
        }
    }
}

void eval_q_inner_loop()
{
    int i;
    int k, r, nbh;
    double x_i, y_i, x_k, y_k;
    double delx, dely, dist, weights;
    double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
    double det, temp;
    double one_by_det;
    double sum_delx_delq1, sum_delx_delq2, sum_delx_delq3, sum_delx_delq4;
    double sum_dely_delq1, sum_dely_delq2, sum_dely_delq3, sum_dely_delq4;
    double q1, q2, q3, q4;

    double temp1, temp2;

    for (int i = 1; i <= local_points; i++)
    {
        x_i = point.x[i];
        y_i = point.y[i];

        sum_delx_sqr = 0.;
        sum_dely_sqr = 0.;
        sum_delx_dely = 0.;

        temp1 = 0.;
        temp2 = 0.;

        sum_delx_delq1 = 0.;
        sum_delx_delq2 = 0.;
        sum_delx_delq3 = 0.;
        sum_delx_delq4 = 0.;

        sum_dely_delq1 = 0.;
        sum_dely_delq2 = 0.;
        sum_dely_delq3 = 0.;
        sum_dely_delq4 = 0.;

        q1 = point.q[1][i];
        q2 = point.q[2][i];
        q3 = point.q[3][i];
        q4 = point.q[4][i];

        for (int k = 1; k <= point.nbhs[i]; k++)
        {
            nbh = point.conn[i][k];

            x_k = point.x[nbh];
            y_k = point.y[nbh];

            delx = x_k - x_i;
            dely = y_k - y_i;

            dist = sqrt(delx * delx + dely * dely);
            weights = pow(dist, power);

            sum_delx_sqr = sum_delx_sqr + delx * delx * weights;
            sum_dely_sqr = sum_dely_sqr + dely * dely * weights;

            sum_delx_dely = sum_delx_dely + delx * dely * weights;

            temp1 = q1 - 0.5 * (delx * point.dq[1][1][i] + dely * point.dq[2][1][i]);
            temp2 = point.q[1][nbh] - 0.5 * (delx * point.dq[1][1][nbh] + dely * point.dq[2][1][nbh]);
            sum_delx_delq1 = sum_delx_delq1 + (weights * delx * (temp2 - temp1));
            sum_dely_delq1 = sum_dely_delq1 + (weights * dely * (temp2 - temp1));

            temp1 = q2 - 0.5 * (delx * point.dq[1][2][i] + dely * point.dq[2][2][i]);
            temp2 = point.q[2][nbh] - 0.5 * (delx * point.dq[1][2][nbh] + dely * point.dq[2][2][nbh]);
            sum_delx_delq2 = sum_delx_delq2 + (weights * delx * (temp2 - temp1));
            sum_dely_delq2 = sum_dely_delq2 + (weights * dely * (temp2 - temp1));

            temp1 = q3 - 0.5 * (delx * point.dq[1][3][i] + dely * point.dq[2][3][i]);
            temp2 = point.q[3][nbh] - 0.5 * (delx * point.dq[1][3][nbh] + dely * point.dq[2][3][nbh]);
            sum_delx_delq3 = sum_delx_delq3 + (weights * delx * (temp2 - temp1));
            sum_dely_delq3 = sum_dely_delq3 + (weights * dely * (temp2 - temp1));

            temp1 = q4 - 0.5 * (delx * point.dq[1][4][i] + dely * point.dq[2][4][i]);
            temp2 = point.q[4][nbh] - 0.5 * (delx * point.dq[1][4][nbh] + dely * point.dq[2][4][nbh]);
            sum_delx_delq4 = sum_delx_delq4 + (weights * delx * (temp2 - temp1));
            sum_dely_delq4 = sum_dely_delq4 + (weights * dely * (temp2 - temp1));
        }

        det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely * sum_delx_dely;
        one_by_det = 1.0 / det;

        point.temp[1][1][i] = one_by_det * (sum_delx_delq1 * sum_dely_sqr - sum_dely_delq1 * sum_delx_dely);
        point.temp[1][2][i] = one_by_det * (sum_delx_delq2 * sum_dely_sqr - sum_dely_delq2 * sum_delx_dely);
        point.temp[1][3][i] = one_by_det * (sum_delx_delq3 * sum_dely_sqr - sum_dely_delq3 * sum_delx_dely);
        point.temp[1][4][i] = one_by_det * (sum_delx_delq4 * sum_dely_sqr - sum_dely_delq4 * sum_delx_dely);
        point.temp[2][1][i] = one_by_det * (sum_dely_delq1 * sum_delx_sqr - sum_delx_delq1 * sum_delx_dely);
        point.temp[2][2][i] = one_by_det * (sum_dely_delq2 * sum_delx_sqr - sum_delx_delq2 * sum_delx_dely);
        point.temp[2][3][i] = one_by_det * (sum_dely_delq3 * sum_delx_sqr - sum_delx_delq3 * sum_delx_dely);
        point.temp[2][4][i] = one_by_det * (sum_dely_delq4 * sum_delx_sqr - sum_delx_delq4 * sum_delx_dely);
    }
}

void eval_update_innerloop()
{
    for (int i = 1; i <= local_points; i++)
    {
        point.dq[1][1][i] = point.temp[1][1][i];
        point.dq[1][2][i] = point.temp[1][2][i];
        point.dq[1][3][i] = point.temp[1][3][i];
        point.dq[1][4][i] = point.temp[1][4][i];
        point.dq[2][1][i] = point.temp[2][1][i];
        point.dq[2][2][i] = point.temp[2][2][i];
        point.dq[2][3][i] = point.temp[2][3][i];
        point.dq[2][4][i] = point.temp[2][4][i];
    }
}

void qtilde_to_primitive(std::vector<double> &qtilde, double &u1, double &u2, double &rho, double &pr)

{
    double beta, temp, temp1, temp2;
    double q1, q2, q3, q4;

    q1 = qtilde[1];
    q2 = qtilde[2];
    q3 = qtilde[3];
    q4 = qtilde[4];

    beta = -q4 * 0.5;

    temp = 0.5 / beta;

    u1 = q2 * temp;
    u2 = q3 * temp;

    temp1 = q1 + beta * (u1 * u1 + u2 * u2);
    temp2 = temp1 - (log(beta) / (gamma_new - 1));

    rho = exp(temp2);
    pr = rho * temp;
}