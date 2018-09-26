#!/usr/bin/bash
#
# Specify name of job allocation
#SBATCH --job-name=cluster_call
#
# Instruct Slurm to connect the batch script's standard output directly
# to the file name specified in the "filename pattern". By default both
# standard output and standard error are directed to the same file
#SBATCH --output=.slurm_%j.out
#
# Specify number of tasks, in general --ntasks should always be equal to 1 (unless you use MPI)
#SBATCH --ntasks=1
#
# Specify CPUs per task, if your application requires n threads, then set --cpus-per-tasks=n (unless MPI is being used)
#SBATCH --cpus-per-task=1
#
# Set a limit on the total run time of the job allocation
#SBATCH --time=1-00:00:00
#
# Specify the memory required per allocated CPU
#SBATCH --mem=20000




Rscript ../exec/clustmut.R -i data/ \
                --glob "*.tsv" \
                --recursive \
                --mode distance \
                -o mela_au_fdr \
                -N 1 \
                -Vvu 
