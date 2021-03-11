#pragma once

#include "data_structure_mod.h"
#include "flux_residual_mod.h"
#include "state_update_mod.h"
#include "q_variables_mod.h"
#include "objective_function_mod.h"
#include "post_processing_mod.h"
#include <iostream>
#include<iomanip>
#include"dump.h"
using namespace std;

void fpi_solver(int temp1)

{

        int i, rk;
        // printf("THis is temp1:%d\n", temp1);
        for (i = 1; i <= max_points; i++)
        {
                for (int j = 1; j <= 4; j++)
                {
                        point.prim_old[j][i] = point.prim[j][i];
                }
        }

        func_delta();
        // printf("THis is temp1:%d\n", temp1);
        //    Perform 4-stage, 3-order SSPRK update
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
        cout<<setprecision(13)<<scientific;
        dump();
        // cout<<"\nq"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {
        //         cout << point.q[i][160] <<  " ";
        // }
        // cout<<"\ndq1"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.dq[1][i][2]<<  " " ;
        // }
        // cout<<"\ndq2"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.dq[2][i][160] <<  " ";
        // }
        // cout<<"\nprim"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.prim[i][1] <<  " ";
        // }
        // cout<<"\nflux_res"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.flux_res[i][1] <<  " ";
        // }
        // cout<<"\nqm1"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.qm[1][i][1]<<  " " ;
        // }
        // cout<<"\nqm2"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.qm[2][i][1] <<  " ";
        // }
        // cout<<"\nprim"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.prim[i][1] <<  " ";
        // }
        // cout<<endl;
        // objective_function();

        res_new = sqrt(sum_res_sqr) / max_points;
        // printf("THis is temp1:%d\n", temp1);
        if (temp1 <= 2 && restart == 0)
        {
                res_old = res_new;
                residue = 0;
        }
        else
        {
                residue = log10(res_new / res_old);
        }
        // cout<<sum_res_sqr<<" "<<res_new<<" "<<res_old<<endl;
        //  Print primal output
        printf("residue(fpi_solver) : %lf\n", residue);
        // if (it % nsave == 0)
        // {
        //         printf("\n");
        //         printf(".............-Saving solution-.............\n");
        //         printf("\n");
        //         print_primal_output();
        // }
}