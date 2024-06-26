#!/bin/bash
#SBATCH --job-name=random_sample
#              d-hh:mm:ss
#SBATCH --time=30-00:00:00
#SBATCH --output=/home/franzf/OpenQuantumComputing_private/OpenQuantumComputing/openquantumcomputing/%j.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32

source activate qiskit

python random_sample_H_overall.py