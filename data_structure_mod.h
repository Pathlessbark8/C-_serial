#pragma once
#include <vector>
#include "parameter_mod.h"
int max_points, local_points, ghost_points;
int wall_points, interior_points, outer_points, shape_points;

//    ghost global indices
int *pghost;

struct points
{

    std::vector<int> original_id;
    std::vector<double> x, y;
    std::vector<int> left, right;
    std::vector<int> flag_1; // stores location of point
    std::vector<int> flag_2; //stores shape point belongs to
    std::vector<int> qtdepth;
    std::vector<double> nx, ny;
    std::vector<int> nbhs;
    std::vector<std::vector<int>> conn; //2-D array

    std::vector<double> min_dist;

    //2-D array
    std::vector<std::vector<double>> prim;
    std::vector<std::vector<double>> prim_old;
    std::vector<std::vector<double>> flux_res;
    std::vector<std::vector<double>> q;
    std::vector<std::vector<double>> U;

    //3-D array
    std::vector<std::vector<std::vector<double>>> dq;
    std::vector<std::vector<std::vector<double>>> qm;
    std::vector<std::vector<std::vector<double>>> temp;

    std::vector<double> entropy, vorticity, vorticity_sqr;
    std::vector<int> xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs;

    std::vector<std::vector<int>> xpos_conn, xneg_conn;
    std::vector<std::vector<int>> ypos_conn, yneg_conn;

    std::vector<double> delta;

    //  Implicit data
    std::vector<std::vector<double>> U_old;
};

points point;

std::vector<int> wall_points_index;
std::vector<int> outer_points_index;
std::vector<int> interior_points_index;
std::vector<int> shape_points_index;

// iterations
int it, itr;

// Flag for time stepping
int rks;
double euler;
double total_loss_stagpressure;
double res_old, res_new, residue, max_res;
double gsum_res_sqr, sum_res_sqr;
int max_res_point;
std::vector<double> Cl, Cd, Cm, cfv, ClCd, vector_cost_func;
double total_entropy, total_enstrophy;
int plen;
int format;

// The parameter CFL is the CFL number for stability ..
double CFL;

int max_iters;

// Unsteady variables
double t, tfinal, dtg;
int timestep;

// Run option: petsc or normal
int runop;

//    The parameter power is used to specify the weights
//    in the LS formula for the derivatives ..
//    power = 0.0d0, -2.0d0, -4.0d0, -6.0d0 ..
//    For example, power = -2.0 implies that
//    power = -2.0 => weights = 1/d^2
//    power = -4.0 => weights = 1/d^4

double power;

//    limiter_flag = 1 => venkatakrishnan limiter
//    limiter_flag = 2 => min-max limiter

int limiter_flag;
double VL_CONST; // Venkatakrishnan limiter constant ..

int restart;

//    Interior points normal flag ..
//    If flag is zero => nx = 0.0 and ny = 1.0

int interior_points_normal_flag;

//   Restart solution parameter
int solution_restart;

//    solution save parameter
int nsave;

//    First order flag
double fo_flag;

//    Objective function
double Cl_flag, Cd_flag, Cm_flag, Cl_Cd_flag, ent_flag, ens_flag;
int obj_flag;

int inner_iterations = 0;

//    No of shapes
int shapes = 1;

void allocate_soln()
{
    point.prim = std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1));
    point.prim_old = std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1));
    point.flux_res = std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1));
    point.U_old = std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1));
    point.q = std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1));
    point.U = std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1));
    point.dq = std::vector<std::vector<std::vector<double>>>(2 + 1, std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1)));
    point.qm = std::vector<std::vector<std::vector<double>>>(2 + 1, std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1)));
    point.temp = std::vector<std::vector<std::vector<double>>>(3 + 1, std::vector<std::vector<double>>(4 + 1, std::vector<double>(max_points + 1)));
    point.entropy = std::vector<double>(max_points + 1);
    point.vorticity = std::vector<double>(max_points + 1);
    point.vorticity_sqr = std::vector<double>(max_points + 1);
    point.xpos_nbhs = std::vector<int>(max_points + 1);
    point.xneg_nbhs = std::vector<int>(max_points + 1);
    point.ypos_nbhs = std::vector<int>(max_points + 1);
    point.yneg_nbhs = std::vector<int>(max_points + 1);
    point.xpos_conn = std::vector<std::vector<int>>(max_points + 1, std::vector<int>(20 + 1));
    point.xneg_conn = std::vector<std::vector<int>>(max_points + 1, std::vector<int>(20 + 1));
    point.ypos_conn = std::vector<std::vector<int>>(max_points + 1, std::vector<int>(20 + 1));
    point.yneg_conn = std::vector<std::vector<int>>(max_points + 1, std::vector<int>(20 + 1));
    point.delta = std::vector<double>(max_points + 1);
    // Cl = std::vector<double>(shapes+1);
    // Cd = std::vector<double>(shapes+1);
    // Cm = std::vector<double>(shapes+1);
    // cfv = std::vector<double>(shapes+1);
    // ClCd = std::vector<double>(shapes+1);
    // vector_cost_func = std::vector<double>(shapes+1);
}
