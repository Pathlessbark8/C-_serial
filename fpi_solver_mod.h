#pragma once

#include "data_structure_mod.h"
#include "flux_residual_mod.h"
#include "state_update_mod.h"
#include "q_variables_mod.h"
#include "objective_function_mod.h"
#include "post_processing_mod.h"
#include <iostream>
#include <iomanip>
#include "dump.h"
using namespace std;

void fpi_solver(int temp1)

{

        int i, rk;
        for (i = 1; i <= max_points; i++)
        {
                for (int j = 0; j < 4; j++)
                {
                        point.prim_old[j][i] = point.prim[j][i];
                }
        }
        // cout<<setprecision(14)<<std::scientific;
        // for(int r=1;r<=4;r++)
        // {
        //         cout<<point.prim[r][0]<<" ";
        // }
        // cout<<endl;
        func_delta();
        // std::cout << "delta " << point.delta[0] << endl;

        // printf("THis is temp1:%d\n", temp1);
        //    Perform 4-stage, 3-order SSPRK update
        for (rk = 1; rk <= rks; rk++)
        {
                // std::cout << "rk : " << rk;
                // std::cout << "\nq" << endl;
                // for (int i = 1; i <= 4; i++)
                // {
                //         std::cout << point.q[i][0] << " ";
                // }
                // std::cout << endl
                //           << endl;

                eval_q_variables();

                eval_q_derivatives();

                for (i = 1; i <= inner_iterations; i++)
                {
                        eval_q_inner_loop();

                        eval_update_innerloop();
                }

                cal_flux_residual();

                state_update(rk);
                // cout << "flux_res" << endl;
                // for (int r = 1; r <= 4; r++)
                // {
                //         cout << point.flux_res[r][0] << " ";
                // }
                // cout << endl;
                // cout<<"Euler :"<<euler<<endl;
                // cout << "\nprim" << endl;
                // for (int i = 1; i <= 4; i++)
                // {

                //         cout << point.prim[i][0] << " ";
                // }
                if (rk == 1)
                {
                        // for(int i=1;i<=4;i++)
                        // {
                        //         std::cout<<point.prim[i][0]<<" ";
                        // }
                        // cout<<endl;
                }
        }

        std::cout << setprecision(13) << scientific;
        // dump();
        // cout<<"rk : "<<rk<<endl;
        // cout<<"\nq"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {
        //         cout << point.q[i][0] <<  " ";
        // }
        // cout<<"\ndq1"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.dq[0][i][1]<<  " " ;
        // }
        // cout<<"\ndq2"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.dq[1][i][160] <<  " ";
        // }
        // cout<<"\nprim"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.prim[i][0] <<  " ";
        // }
        // cout<<"\nflux_res"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.flux_res[i][0] <<  " ";
        // }
        // cout<<"\nqm1"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.qm[0][i][0]<<  " " ;
        // }
        // cout<<"\nqm2"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.qm[1][i][0] <<  " ";
        // }
        // cout<<"\nprim"<<endl;
        // for (int i = 1; i <= 4; i++)
        // {

        //         cout << point.prim[i][0] <<  " ";
        // }
        // cout<<endl;
        // objective_function();

        res_new = sqrt(sum_res_sqr) / max_points;
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
        cout << "Iterations "<< it << " Residue : " << residue << "\n";//" res_new :" << res_new << " res_old :" << res_old << endl;
}