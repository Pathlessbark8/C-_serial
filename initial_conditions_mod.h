#pragma once
    #include"data_structure_mod.h"
    #include"parameter_mod.h"


    void initial_conditions()

        {
        int k,i;

        if(restart == 0) 

            {setup_case_parameters();

            for(k=1;k<= max_points;k++)
                {point.prim[1][k] = q_init[1];
                point.prim[2][k] = q_init[2];
                point.prim[3][k] = q_init[3];
                point.prim[4][k] = q_init[4];
            }
            }
        
        else if(restart == 1) 
    
           { // if(rank == 0) then
            //     write(*,*)'............-Using restart file-...........'
            //     write(*,*)
            // end if

             setup_case_parameters();
            
             //restart_sol();

            //  update_begin_prim_ghost();
            //  update_end_prim_ghost();
        
      }

}

//     void restart_sol()
 
//       {

//         int i, dummy;
//         char sfile[64],itos[64];

//         if(proc==1) 
//             sfile = "restart/sol.dat"
//         else
//             sfile = 'restart/'//'sol-'//trim(itos(4,rank))//'.dat'
//         end if

//         OPEN(UNIT=515,FILE=trim(sfile),form='formatted', action="read")

//         read(515,*)dummy, itr, res_old

//         do i = 1, local_points
//             read(515,*)dummy, dummy, dummy, dummy,&
//                 dummy, point.prim(1,i), point.prim(2,i), point.prim(3,i),&
//                 point.prim(4,i)
//         end do

//         close(515)

// }
