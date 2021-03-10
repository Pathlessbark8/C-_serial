        #include"parameter_mod.h"
        #include"data_structure_mod.h"
        #include"q_lskum_mod.h"
        #include"compute_force_coeffs_mod.h"


void main(){
        double totaltime,runtime;

        if(ierr /= 0)
        {
        printf("Unable to initialize PETSc");
            exit(0);    
        } 
        // stop "Unable to initialize PETSc"
         MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr)
         MPI_Comm_size(PETSC_COMM_WORLD, proc, ierr)
        if(rank==0) then
                 execute_command_line('mkdir -p solution')
                 execute_command_line('mkdir -p cp')
        end if
        
        totaltime = MPI_Wtime()

        if(rank == 0) then
                write(*,*)
                write(*,*)'%%%%%%%%%%%%%-MPI Meshfree Code-%%%%%%%%%%%'
                write(*,*)
        end if

!       Read the case file

        if(rank == 0) then
                write(*,*)'%%%%%%%%%%%-Reading the case file-%%%%%%%%%'
                write(*,*)
        end if
        
         readcase()

!	Reading the input data ..

        if(rank == 0) then
                write(*,*)
                write(*,*)'%%%%%%%%%%%%-Reading HDF5 point data-%%%%%%%%%%%'
                write(*,*)
        end if
        
         read_hdf5input_point_data()

        !  read_input_point_data()
        
!       Allocate solution variables

         allocate_soln()

!       Initialize Petsc vectors

        if(proc .ne. 1) init_petsc()
        if(proc == 1) plen = max_points
        if(rank == 0) then
                write(*,*) 'Number of points:         ', plen
                write(*,*)
        end if
        
!	Assign the initial conditions for the primitive variables ..	
         initial_conditions()
        if(rank == 0) then
                write(*,*)'%%%%%%%%%%%-Solution initialised-%%%%%%%%%%'
                write(*,*)
        end if

!	Primal fixed point iterative solver ..
        runtime = MPI_Wtime()
        if(runop == 1)then
                if(rank == 0) then
                        write(*,*)'%%%%%%%%%-Using inbuilt solvers-%%%%%%%%%%%'
                        write(*,*)
                end if
                 q_lskum()
        end if
        runtime = MPI_Wtime() - runtime

!       Save solution one last time

        !  print_primal_output()

!       destroy petsc vectors and deallocate point/solution vectors
         dest_petsc()
         deallocate_soln()
         dealloc_points()

        totaltime = MPI_Wtime() - totaltime

        if(rank == 0) then
                write(*,*)
                write(*,*) '%%%%%%%%%%%-Simulation finished-%%%%%%%%%%%'
                write(*,*) 'Run time:  ',runtime,'seconds'
                write(*,*) 'Total time:',totaltime,'seconds'
        end if

!       stop petsc
         PetscFinalize(ierr)
}