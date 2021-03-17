#!/usr/bin/env bash
#SBATCH --time=5:30:00
#SBATCH --partition=hsw_v100_32g
#SBATCH --job-name="base_128"
spack load nvhpc%gcc@10.2.0
spack load gcc@9.3.0
rm partGrid
ln -s /home/aniln/grids/c++/40M partGrid
nvcc --compiler-options -mcmodel=medium -std=c++14 main.cu
CUDA_VISIBLE_DEVICES=0 ./a.out
#CUDA_VISIBLE_DEVICES=0 ncu -o 128block_20M_mill_final --set full ./a.out

