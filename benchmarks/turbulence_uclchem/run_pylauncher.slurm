#!/bin/bash
#SBATCH -N 1                      # Request 2 nodes (56 cores each = 112 total)
#SBATCH -p gg
#SBATCH -t 01:00:00
#SBATCH -A AST24021
#SBATCH -J uclchem_turbulence_benchmark
#SBATCH -o /work/09338/carlos9/vista/slurmjobinfo/job.o%j

module load gcc/13.2.0
module load pylauncher

source /work/09338/carlos9/vista/venv/bin/activate

python /work/09338/carlos9/vista/UCLCHEM/benchmarks/turbulence/mylauncher.py
