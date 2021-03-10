#pragma once
        
       #include"data_structure_mod.h"
       #include"flux_residual_mod.h"
        #include"state_update_mod.h"
        #include"q_variables_mod.h"
        #include"objective_function_mod.h"
        #include"post_processing_mod.h" 

        void fpi_solver(int t)

{
                
                int i, rk;

                for(i =1;i<= max_points;i++)
                {
                    for(int j=1;j<=4;j++)
                    {
                        point.prim_old[j][i] = point.prim[j][i];
                }
                }

                 func_delta();

            //    Perform 4-stage, 3-order SSPRK update
                for( rk = 1;rk<= rks;rk++)
                {
                         eval_q_variables();
                
                         eval_q_derivatives();

                                for( i = 1;i<= inner_iterations;i++)
{
                                         eval_q_inner_loop();

                                         eval_update_innerloop();
                                        
}
                
                         cal_flux_residual();

                         state_update(rk);

}
                
                 objective_function();

                res_new = sqrt(sum_res_sqr)/max_points;

                if(t < 2 && restart == 0)
                        {res_old = res_new;
                        residue = 0;}
                else 
                       { residue = log10(res_new/res_old);}

                //  Print primal output
                if(it%nsave==0) {
                        printf("\n");
                        printf(".............-Saving solution-.............\n");
                        printf("\n");
                         print_primal_output();
                }

  }