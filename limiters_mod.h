#pragma once

#include "data_structure_mod.h"
#include <vector>

// The following two subroutines are used for min-max limiter ..

void max_q_value(int i, double maxi[])
{
    int j, k, r;

    for (int j = 0; j < 4; j++)
    {
        maxi[j] = point.q[j][i];
    }

    for (j = 1; j <= point.nbhs[i]; j++)
    {
        k = point.conn[i][j];
        for (int r = 0; r < 4; r++)
        {
            if (maxi[r] < point.q[r][k])
            {
                maxi[r] = point.q[r][k];
            }
        }
    }
}

void min_q_value(int i, double mini[])
{

    int j, k, r;

    for (int j = 0; j <4; j++)
    {
        mini[j] = point.q[j][i];
    }

    for (j = 1; j <= point.nbhs[i]; j++)
    {
        k = point.conn[i][j];
        for (int r = 0; r < 4; r++)
        {
            if (mini[r] > point.q[r][k])
            {
                mini[r] = point.q[r][k];
            }
        }
    }
}

// The following subroutines are used for venkatakrishnan limiter ..

void venkat_limiter(double qtilde[], double phi[], int k)
{

    int r;
    double q, del_neg, del_pos;
    double max_q, min_q, ds, epsi, num, den, temp;

    for (r = 0; r < 4; r++)
    {
        q = point.q[r][k];
        del_neg = qtilde[r] - q;
        if (abs(del_neg) <= 10e-6)
        {
            phi[r] = 1.0;
        }
        else if (abs(del_neg) > 10e-6)
        {
            if (del_neg > 0)
            {
                del_pos = point.qm[0][r][k] - q;
            }

            else if (del_neg < 0)
            {
                del_pos = point.qm[1][r][k] - q;
            }

            epsi = VL_CONST * point.min_dist[k];
            epsi = pow(epsi, 3);

            num = (del_pos * del_pos) + (epsi * epsi); // Numerator ..
            num = num * del_neg + 2.0 * del_neg * del_neg * del_pos;

            den = del_pos * del_pos + 2.0 * del_neg * del_neg; // Denominator ..
            den = den + del_neg * del_pos + epsi * epsi;
            den = den * del_neg;

            temp = num / den;

            if (temp < 1)
            {
                phi[r] = temp;
            }
            else
            {
                phi[r] = 1;
            }
        }
    }
}
