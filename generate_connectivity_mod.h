#pragma once

#include "data_structure_mod.h"

int return_index(int pt_to_search, int host_index){
    unsigned int n = sizeof(point.conn[0])/sizeof(point.conn[0][0]);
    for (int i = 1; i < n; i++){
        if (point.conn[host_index][i] == pt_to_search){
            return i;
        }
    }
    exit(2);
}

void get_interior_neighbours(int i, double nx, double ny)

{

    double xi, yi, xk, yk;
    double delx, dely, dels, deln;
    double tx, ty;
    int r, count, nbh;

    xi = point.x[i];
    yi = point.y[i];

    tx = ny;
    ty = -nx;

    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;

    for (r = 1; r <= point.nbhs[i]; r++)
    {
        nbh = point.conn[i][r];
        xk = point.x[nbh];
        yk = point.y[nbh];

        delx = xk - xi;
        dely = yk - yi;

        dels = delx * tx + dely * ty;
        deln = delx * nx + dely * ny;

        if (dels <= 0.0)

        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;

            count = point.xpos_nbhs[i];
            point.xpos_conn[i][count] = return_index(nbh, i);
        }

        if (dels >= 0.0)

        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;

            count = point.xneg_nbhs[i];
            point.xneg_conn[i][count] = return_index(nbh, i);
        }

        if (deln <= 0.0)

        {
            point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;

            count = point.ypos_nbhs[i];
            point.ypos_conn[i][count] = return_index(nbh, i);
        }

        if (deln >= 0.0)

        {
            point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;

            count = point.yneg_nbhs[i];
            point.yneg_conn[i][count] = return_index(nbh, i);
        }
    }

    // if(point.xpos_nbhs[i] == 0)
    //     print*,"WARNING!!! xpos zero for interior point number:", i,", rank:", rank
    // else if(point.xneg_nbhs[i] == 0)
    //     print*,"WARNING!!! xneg zero for interior point number:", i,", rank:", rank
    // else if(point.ypos_nbhs[i] == 0)
    //     print*,"WARNING!!! ypos zero for interior point number:", i,", rank:", rank
    // else if(point.yneg_nbhs[i] == 0)
    //     print*,"WARNING!!! yneg zero for interior point number:", i,", rank:", rank
    // end if
}
void check_condition_number(int i, double nx, double ny)
{
    // Use lapack_example_aux, Only: nagf_file_print_matrix_real_gen
    // Use lapack_interfaces, Only: dbdsqr, dgebrd, dlacpy, dorgbr
    // Use lapack_precision, Only: dp
}

void get_wall_boundary_neighbours(int i, double nx, double ny)
{

    double xi, yi, xk, yk;
    double delx, dely, dels, deln;
    double tx, ty;
    int r, count, nbh;

    xi = point.x[i];
    yi = point.y[i];

    tx = ny;
    ty = -nx;

    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.yneg_nbhs[i] = 0;

    for (r = 1; r <= point.nbhs[i]; r++)
    {
        nbh = point.conn[i][r];

        xk = point.x[nbh];
        yk = point.y[nbh];

        delx = xk - xi;
        dely = yk - yi;

        dels = delx * tx + dely * ty;
        deln = delx * nx + dely * ny;

        if (dels <= 0.0)

        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;

            count = point.xpos_nbhs[i];
            point.xpos_conn[i][count] = return_index(nbh, i);
        }

        if (dels >= 0.0)

        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;

            count = point.xneg_nbhs[i];
            point.xneg_conn[i][count] = return_index(nbh, i);
        }

        point.yneg_nbhs[i] = point.yneg_nbhs[i] + 1;

        count = point.yneg_nbhs[i];
        point.yneg_conn[i][count] = return_index(nbh, i);
    }

    // if(point.xpos_nbhs[i] == 0)
    //     print*,"WARNING!!! xpos zero for wall point number:", i,", rank:", rank
    // elseif(point.xneg_nbhs[i] == 0)
    //     print*,"WARNING!!! xneg zero for wall point number:", i,", rank:", rank
    // elseif(point.yneg_nbhs[i] == 0)
    //     print*,"WARNING!!! yneg zero for wall point number:", i,", rank:", rank
    // end if
}

void get_outer_boundary_neighbours(int i, double nx, double ny)

{

    double xi, yi, xk, yk;
    double delx, dely, dels, deln;
    double tx, ty;
    int r, count, nbh;

    xi = point.x[i];
    yi = point.y[i];

    tx = ny;
    ty = -nx;

    point.xpos_nbhs[i] = 0;
    point.xneg_nbhs[i] = 0;
    point.ypos_nbhs[i] = 0;

    for (r = 1; r <= point.nbhs[i]; r++)

    {
        nbh = point.conn[i][r];

        xk = point.x[nbh];
        yk = point.y[nbh];

        delx = xk - xi;
        dely = yk - yi;

        dels = delx * tx + dely * ty;
        deln = delx * nx + dely * ny;

        if (dels <= 0.0)

        {
            point.xpos_nbhs[i] = point.xpos_nbhs[i] + 1;

            count = point.xpos_nbhs[i];
            point.xpos_conn[i][count] = return_index(nbh, i);
        }

        if (dels >= 0.0)

        {
            point.xneg_nbhs[i] = point.xneg_nbhs[i] + 1;

            count = point.xneg_nbhs[i];
            point.xneg_conn[i][count] = return_index(nbh, i);
        }

        point.ypos_nbhs[i] = point.ypos_nbhs[i] + 1;

        count = point.ypos_nbhs[i];
        point.ypos_conn[i][count] = return_index(nbh, i);
    }

    // if(point.xpos_nbhs[i] == 0)
    //     print*,"WARNING!!! xpos zero for outer point number:", i,", rank:", rank
    // elseif(point.xneg_nbhs[i] == 0)
    //     print*,"WARNING!!! xneg zero for outer point number:", i,", rank:", rank
    // elseif(point.ypos_nbhs[i] == 0)
    //     print*,"WARNING!!! ypos zero for outer point number:", i,", rank:", rank
    // end if
}

void generate_connectivity()

{

    double nx, ny;
    for (int i = 1; i <= max_points; i++)
    {
        if (point.flag_1[i] == 0)
        {
            nx = point.nx[i];
            ny = point.ny[i];
            get_wall_boundary_neighbours(i, nx, ny);
        }
        else if (point.flag_1[i] == 1)
        {
            nx = point.nx[i];
            ny = point.ny[i];
            get_interior_neighbours(i, nx, ny);
            check_condition_number(i, nx, ny);
        }
        else
        {
            nx = point.nx[i];
            ny = point.ny[i];
            get_outer_boundary_neighbours(i, nx, ny);
        }
    }
}
