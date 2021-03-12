#pragma once
       
    #include "data_structure_mod.h"

        void compute_enstrophy()
 {           
            int i, k, r, nbh;
            double  x_i, y_i, x_k, y_k;
			double  delx, dely, dist, weights;
			double sum_delx_sqr, sum_dely_sqr, sum_delx_dely;
			double sum_delx_delu1, sum_delx_delu2, sum_dely_delu1, sum_dely_delu2;
			double  det;
			double  one_by_det;
			double  du1_dy, du2_dx, temp;
            double  gtotal_enstrophy;

            total_enstrophy = 0;

            for(i = 1;i<= local_points;i++){
            

                x_i = point.x[i];
                y_i = point.y[i];

                sum_delx_sqr = 0.;
                sum_dely_sqr = 0.;
                sum_delx_dely = 0.;

                sum_delx_delu1 = 0.;
                sum_dely_delu1 = 0.;
                sum_delx_delu2 = 0.;
                sum_dely_delu2 = 0.;

                for(k = 1;k<= point.nbhs[i];k++){

                    nbh = point.conn[i][k];

                    x_k = point.x[nbh];
                    y_k = point.y[nbh];

                    delx = x_k - x_i;
                    dely = y_k - y_i;

                    dist = sqrt(delx*delx + dely*dely);
                    weights = pow(dist,power);

                    sum_delx_sqr = sum_delx_sqr + delx*delx*weights;
                    sum_dely_sqr = sum_dely_sqr + dely*dely*weights;

                    sum_delx_dely = sum_delx_dely + delx*dely*weights;

                    sum_delx_delu1 = sum_delx_delu1 + weights*delx*(point.prim[1][nbh] - point.prim[1][i]);
                    sum_delx_delu2 = sum_delx_delu2 + weights*delx*(point.prim[2][nbh] - point.prim[2][i]);
                    sum_dely_delu1 = sum_dely_delu1 + weights*dely*(point.prim[1][nbh] - point.prim[1][i]);
                    sum_dely_delu2 = sum_dely_delu2 + weights*dely*(point.prim[2][nbh] - point.prim[2][i]);

}

                det = sum_delx_sqr*sum_dely_sqr - sum_delx_dely*sum_delx_dely;
                one_by_det = 1.0/det;

                du2_dx = (sum_delx_delu2*sum_dely_sqr - sum_dely_delu2*sum_delx_dely)*one_by_det;
                du1_dy = (sum_dely_delu1*sum_delx_sqr - sum_delx_delu1*sum_delx_dely)*one_by_det;

                temp = du2_dx - du1_dy;

                point.vorticity[i] = temp;

                point.vorticity_sqr[i] = temp*temp;

                total_enstrophy = total_enstrophy + point.vorticity_sqr[i];

            }

}

