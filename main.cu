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
    fin.open("partGrid", ios::in);
    fin >> local_points;
    for (int i = 1; i <= max_points; i++)
    {
        fin >> point.x[i] >> point.y[i] >> point.left[i] >> point.right[i] >> point.flag_1[i] >> point.flag_2[i] >> point.min_dist[i] >> point.nbhs[i];
        for (int r = 1; r <= point.nbhs[i]; r++)
        {
            fin >> point.conn[i][r];
        }
    }
    allocate_soln();
    local_points = max_points;
    compute_normals();
    generate_connectivity();
    initial_conditions();
    q_lskum();
}