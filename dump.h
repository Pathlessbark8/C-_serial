#pragma once

#include "data_structure_mod.h"
#include <fstream>
#include <iomanip>
using namespace std;
void dump()
{
    fstream fout1, fout2, fout3, fout4,fout5;
    fout1.open("dump_q", ios::out);
    fout2.open("dump_dq1", ios::out);
    fout3.open("dump_flux_res", ios::out);
    fout4.open("dump_dq2", ios::out);
    fout5.open("dump_prim", ios::out);
    fout1 << setprecision(14) << std::scientific;
    fout2 << setprecision(14) << std::scientific;
    fout3 << setprecision(14) << std::scientific;
    fout4 << setprecision(14) << std::scientific;
    fout5 << setprecision(14) << std::scientific;
    for (int i = 1; i <= max_points; i++)
    {
        fout1 << i << "   " ;
        for (int j = 1; j <= 4; j++)
        {
            fout1  << point.q[j][i] << " ";
        }
        fout1 << endl;
        fout2 << i << "   " ;
        for (int j = 1; j <= 4; j++)
        {
            fout2  << point.dq[0][j][i] << " ";
        }
        fout2 << endl;
        fout3 << i << "   " ;
        for (int j = 1; j <= 4; j++)
        {
            fout3<< point.flux_res[j][i] << " ";
        }
        fout3 << endl;
        fout4 << i << "   " ;
        for (int j = 1; j <= 4; j++)
        {
            fout4  << point.dq[1][j][i] << " ";
        }
        fout4 << endl;
        fout5 << i << "   " ;
        for (int j = 1; j <= 4; j++)
        {
            fout5  << point.prim[j][i] << " ";
        }
        fout5 << endl;
    }
    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
    fout5.close();
}