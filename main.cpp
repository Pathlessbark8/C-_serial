#include "split_fluxes_mod.h"
#include "quadrant_fluxes_mod.h"
#include "iostream"
#include <fstream>
#include "data_structure_mod.h"
#include <iomanip>
#include "initial_conditions_mod.h"
#include "q_lskum_mod.h"
#include "post_processing_mod.h"

using namespace std;

int main()
{

    fstream fin;
    fin.open("/opt/grids/legacy/part/partGrid9600", ios::in);
    fin >> max_points;
    point.x = std::vector<double>(max_points + 1);
    point.y = std::vector<double>(max_points + 1);
    point.left = std::vector<int>(max_points + 1);
    point.right = std::vector<int>(max_points + 1);
    point.flag_1 = std::vector<int>(max_points + 1);
    point.flag_2 = std::vector<int>(max_points + 1);
    point.nbhs = std::vector<int>(max_points + 1);
    point.min_dist = std::vector<double>(max_points + 1);
    point.conn = std::vector<std::vector<int>>(max_points + 1, std::vector<int>(20 + 1));
    point.nx = std::vector<double>(max_points + 1);
    point.ny = std::vector<double>(max_points + 1);
    for (int i = 1; i <= max_points; i++)
    {
        fin >> point.x[i] >> point.y[i] >> point.left[i] >> point.right[i] >> point.flag_1[i] >> point.flag_2[i] >> point.min_dist[i] >> point.nbhs[i];
        for (int r = 1; r <= point.nbhs[i]; r++)
        {
            fin >> point.conn[i][r];
        }
        if (point.flag_1[i] == 0)
        {
            wall_points = wall_points + 1;
            wall_points_index.push_back(i);
        }
        else if (point.flag_1[i] == 1)
        {
            interior_points = interior_points + 1;
            interior_points_index.push_back(i);
        }
        else if (point.flag_1[i] == 2)
        {
            outer_points = outer_points + 1;
            outer_points_index.push_back(i);
        }

        if (point.flag_2[i] > 0)
        {
            shape_points = shape_points + 1;
        }
        // cout << point.x[i] << " " << point.y[i] << " " << point.left[i] << " " << point.right[i] << " " << point.flag_1[i] << " " << point.flag_2[i] << " " << point.min_dist[i] << " " << point.nbhs[i] << endl;
    }
    allocate_soln();
    local_points = max_points;
    compute_normals();
    generate_connectivity();
    // cout <<"NBHS\n"<< point.xpos_nbhs[1]<<" "<< point.xneg_nbhs[1]<<" "<< point.ypos_nbhs[1]<<" "<< point.yneg_nbhs[1]<<endl;
    // cout <<"NBHS\n"<< point.nbhs[1]<<endl;
    initial_conditions();
    q_lskum();

    // print_primal_output();
    //    cout<<local_points<<endl;
    // double arr[5];
    // double t_nx, t_ny, t_u1,t_u2, t_rho,t_pr;
    // fstream fin,fout;
    // fin.open("data.dat",ios::in);
    // fout.open("validate.dat",ios::out);
    // while(!fin.eof())
    // {
    //     fin>>arr[1]>>arr[2]>>arr[3]>>arr[4];
    //     fin>>t_nx>> t_ny>> t_u1>>t_u2>> t_rho>>t_pr;
    //     if(!fin.eof())
    //     {
    //         // cout<<arr[1]<<arr[2]<<arr[3]<<arr[3]<<endl;
    //         // cout<<t_nx<< t_ny<< t_u1<<t_u2<< t_rho<<t_pr<<endl;
    //         // for(int i=1;i<5;i++)
    //         // {
    //         //     cout<<arr[i]<<" ";
    //         // }
    //         // cout<<endl;
    //         flux_quad_GxI(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
    //         // flux_Gxn(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
    //         // flux_Gyp(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
    //         // flux_Gyn(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
    //         for(int i=1;i<5;i++)
    //         {
    //             fout<<setw(25)<<setprecision(14)<<scientific <<arr[i]<<" ";
    //         }
    //         fout<<endl;
    //     }
    // }
    // fin.close();
    // fout.close();
}