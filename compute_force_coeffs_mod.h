#pragma once

#include "data_structure_mod.h"
#include <vector>
#include <fstream>

void compute_cl_cd_cm()
{

    int j;
    int l, m, r;
    double cp, temp;
    double lx, ly, mx, my, rx, ry;
    double ds1, ds2, ds;

    // std::vector<double> H(shapes + 1, 0), V(shapes + 1, 0), pitch_mom(shapes + 1, 0);
    double H[shapes+1]={};
    double V[shapes+1]={};
    double pitch_mom[shapes+1]={};
    double nx, ny;
    // char cp_file[] = "cp-file";

    std::fstream fout;
    fout.open("cp_file", std::ios::out);

    temp = 0.5 * rho_inf * mach * mach;

    for (j = 1; j <= shape_points; j++)
    {
        m = j;
        r = point.right[m];
        l = point.left[m];

        lx = point.x[l];
        ly = point.y[l];

        mx = point.x[m];
        my = point.y[m];

        rx = point.x[r];
        ry = point.y[r];

        ds1 = pow((mx - lx), 2) + pow((my - ly), 2);
        ds1 = sqrt(ds1);

        ds2 = pow((rx - mx), 2) + pow((ry - my), 2);
        ds2 = sqrt(ds2);

        ds = 0.5 * (ds1 + ds2);

        nx = point.nx[m];
        ny = point.ny[m];

        cp = point.prim[3][m] - pr_inf;
        cp = -cp / temp;

        fout << point.flag_2[m], point.x[m], cp;

        H[point.flag_2[m]] = H[point.flag_2[m]] + cp * nx * ds;
        V[point.flag_2[m]] = V[point.flag_2[m]] + cp * ny * ds;

        pitch_mom[point.flag_2[m]] = pitch_mom[point.flag_2[m]] + (-cp * ny * ds * (mx - 0.25) + cp * nx * ds * (my));
    }

    for (int j = 0; j < shapes; j++)
    {
        Cl[j] = V[j] * cos(theta) - H[j] * sin(theta);
        Cd[j] = H[j] * cos(theta) + V[j] * sin(theta);
        // Cm = pitch_mom;
        std::copy(std::begin(pitch_mom), std::end(pitch_mom), std::begin(Cm));
    }
    fout.close();
}
