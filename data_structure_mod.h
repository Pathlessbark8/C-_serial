#pragma once
#include <vector>
#include "parameter_mod.h"

#define max_points 625000
#define pi 3.141592653589793238462643383279502884

int local_points, ghost_points;
int wall_points, interior_points, outer_points, shape_points;

//    ghost global indices
int *pghost;

struct points
{
  double x[max_points + 1], y[max_points + 1];
  int left[max_points + 1], right[max_points + 1];
  int flag_1[max_points + 1]; // stores location of point
  int flag_2[max_points + 1]; //stores shape point belongs to
  double nx[max_points + 1], ny[max_points + 1];
  int nbhs[max_points + 1];
  int conn[max_points + 1][15]; //2-D array

  double min_dist[max_points + 1];

  double prim[4][max_points + 1];
  double prim_old[4][max_points + 1];
  double flux_res[4][max_points + 1];
  double q[4][max_points + 1];
  double U[4][max_points + 1];
  double dq[2][4][max_points];
  double qm[2][4][max_points];
  double temp[3][4][max_points];

  int xpos_nbhs[max_points + 1], xneg_nbhs[max_points + 1], ypos_nbhs[max_points + 1], yneg_nbhs[max_points + 1];

  u_short xpos_conn[max_points][15], xneg_conn[max_points][15];
  u_short ypos_conn[max_points][15], yneg_conn[max_points][15];

  double delta[max_points];
  double U_old[4][max_points + 1];
};

points point;

std::vector<int> wall_points_index;
std::vector<int> outer_points_index;
std::vector<int> interior_points_index;
std::vector<int> shape_points_index;

//CUDA Specific

int threads_per_block = 128;

// iterations
int it = 0, itr = 3;

// Flag for time stepping
int rks = 4;
double euler = 1.0;
double total_loss_stagpressure;
double res_old = 0, res_new = 0, residue, max_res = 0;
double sum_res_sqr = 0;
int max_res_point = 0;
//    No of shapes
const int shapes = 1;
double Cl[max_points + 1], Cd[max_points + 1], Cm[max_points + 1], cfv[max_points + 1], ClCd[max_points + 1], vector_cost_func[max_points + 1];
double total_entropy, total_enstrophy;
int plen;
int format;

// The parameter CFL is the CFL number for stability ..
double CFL = 0.00001;

int max_iters = 10;

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

int power = 0;

//    limiter_flag = 1 => venkatakrishnan limiter
//    limiter_flag = 2 => min-max limiter

int limiter_flag;
double VL_CONST = 100; // Venkatakrishnan limiter constant ..

int restart = 0;

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

int inner_iterations = 3;

void allocate_soln()
{

}