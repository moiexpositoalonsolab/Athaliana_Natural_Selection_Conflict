#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=DPB
#SBATCH --job-name=pcs_FT16_Growth_rate
#SBATCH --output=pcs_FT16_Growth_rate.slurm.log
../gemma -bfile ../1001gbi -miss 0.05 -maf 0.05 -r2 1 -k ../1001gbi.sXX.txt -lmm 4 -n 1 10 -c 1001pcs1_5.eigenvec -o mGWAS_pcs_FT16_Growth_rate
