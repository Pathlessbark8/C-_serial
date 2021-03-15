#pragma once

#include "data_structure_mod.h"
#include <fstream>
#include <iostream>

// void print_primal_output()

// {
//         int i;
//         double sos, vel_mag, mach_number;
//         std::fstream fout;
//         fout.open("output.dat", std::ios::out);

//         fout << max_points << " " << it << " " << res_old;
//         for (i = 1; i <= max_points; i++)
//         {
//                 sos = sqrt(gamma_new * point.prim[3][i] / point.prim[0][i]);
//                 vel_mag = point.prim[1][i] * point.prim[1][i] + point.prim[2][i] * point.prim[2][i];

//                 mach_number = sqrt(vel_mag / sos);

//                 fout << point.x[i] << " " << point.y[i] << " " << point.flag_1[i] << " "
//                      << " " << point.prim[0][i] << " " << point.prim[1][i] << " " << point.prim[2][i] << " " << point.prim[3][i] << " " << mach_number << " " << point.entropy[i];
//         }
//         fout.close();
// }