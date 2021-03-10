
#pragma once
#include "data_structure_mod.h"

void compute_normals()

{

    double lx, ly, mx, my, rx, ry;
    double nx1, nx2, ny1, ny2, nx, ny;
    double det;

    int i, j, k, l, m, r;

    // Finding the normals for the points on the shapes ..

    for (i = 1; i <= wall_points; i++)
    {

        m = wall_points_index[i];
        l = point.left[m];
        r = point.right[m];

        lx = point.x[l];
        ly = point.y[l];

        mx = point.x[m];
        my = point.y[m];

        rx = point.x[r];
        ry = point.y[r];

        nx1 = my - ly;
        nx2 = ry - my;

        ny1 = mx - lx;
        ny2 = rx - mx;

        nx = 0.5 * (nx1 + nx2);
        ny = 0.5 * (ny1 + ny2);

        det = sqrt(nx * nx + ny * ny);

        nx = -nx / det;
        ny = ny / det;

        point.nx[m] = nx;
        point.ny[m] = ny;
    }

    // Finding the normals for the outer boundary points ..

    for (i = 1; i <= outer_points; i++)
    {
        m = outer_points_index[i];
        l = point.left[m];
        r = point.right[m];

        lx = point.x[l];
        ly = point.y[l];

        mx = point.x[m];
        my = point.y[m];

        rx = point.x[r];
        ry = point.y[r];

        nx1 = my - ly;
        nx2 = ry - my;

        ny1 = mx - lx;
        ny2 = rx - mx;

        nx = 0.5 * (nx1 + nx2);
        ny = 0.5 * (ny1 + ny2);

        det = sqrt(nx * nx + ny * ny);

        nx = -nx / det;
        ny = ny / det;

        point.nx[m] = nx;
        point.ny[m] = ny;
    }

    if (interior_points_normal_flag == 0 && format != 2)
    {
        for (i = 1; i <= interior_points; i++)
        {
            k = interior_points_index[i];
            point.nx[k] = 0.;
            point.ny[k] = 1.;
        }
    }
    else if (interior_points_normal_flag == 1 && format != 2)
    {
        for (i = 1; i <= interior_points; i++)
        {
            k = interior_points_index[i];
            point.nx[k] = 1.;
            point.ny[k] = 0.;
        }
    }
}
