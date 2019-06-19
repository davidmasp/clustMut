#!/usr/bin/bash

#SBATCH --job-name=TEST_clustmut

#SBATCH --output=.slurm_%j.out

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=1

#SBATCH --time=00:60:00

#SBATCH --mem=10000


srun Rscript ../exec/clustmut_distance.R -i data/ \
                --glob "*TCGA-CG-5723-01A-11D-1600-08_randomized_FN_0.75_simulated_clust_rmdup.tsv_w500000.randomized.tsv" \
                --recursive \
                -b "tests_exec/data/TCGA-CG-5723-01A-11D-1600-08_mutations_strand2.txt" \
                -o mela_au_fdr \
                -N 1 \
                -Vvu
#
# srun Rscript ../exec/clustmut_distance.R -i data/ \
#                 --glob "*.tsv" \
#                 --recursive \
#                 -o mela_au_fdr \
#                 -N 1 \
#                 -Vvu
