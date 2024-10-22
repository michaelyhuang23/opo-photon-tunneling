#!/bin/bash

#SBATCH --cpus-per-task=20
#SBATCH -o sol_compare_output
#SBATCH --job-name=clip

module load julia/1.10.1

julia sol_compare.jl
