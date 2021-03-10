#include "split_fluxes_mod.h"
#include "quadrant_fluxes_mod.h"
#include "iostream"
#include <fstream>

#include<iomanip>

using namespace std;

int main(){

    double arr[5];
    double t_nx, t_ny, t_u1,t_u2, t_rho,t_pr;
    fstream fin,fout;
    fin.open("data.dat",ios::in);
    fout.open("validate.dat",ios::out);
    while(!fin.eof())
    {
        fin>>arr[1]>>arr[2]>>arr[3]>>arr[4];
        fin>>t_nx>> t_ny>> t_u1>>t_u2>> t_rho>>t_pr;
        if(!fin.eof())
        {
            // cout<<arr[1]<<arr[2]<<arr[3]<<arr[3]<<endl;
            // cout<<t_nx<< t_ny<< t_u1<<t_u2<< t_rho<<t_pr<<endl;
            // for(int i=1;i<5;i++)
            // {
            //     cout<<arr[i]<<" ";
            // }
            // cout<<endl;
            flux_quad_GxI(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
            // flux_Gxn(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
            // flux_Gyp(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
            // flux_Gyn(arr,t_nx, t_ny, t_u1,t_u2, t_rho,t_pr);
            for(int i=1;i<5;i++)
            {
                fout<<setw(25)<<setprecision(14)<<scientific <<arr[i]<<" ";
            }
            fout<<endl;
        }
    }
    fin.close();
    fout.close();
}