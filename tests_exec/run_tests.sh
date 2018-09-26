#!/usr/bin/bash

# run tests

# from https://stackoverflow.com/questions/592620
if hash sbatch 2>/dev/null; then
    echo "Running exec test of clustmut in slurm mode"

    sbatch --wait clean.sh

    sbatch run_test_lfdr.sh 
    sbatch --wait run_test_sFDR.sh

    srun --job-name "test_check" Rscript check_output_status.R

else
    echo "Running exec test of clustmut in local mode"

    sh clean.sh
    sh run_test_fdr.sh
    sh run_test_FDR.sh
    
    srun --jobname "test_check" Rscript check_output_status.R
fi


