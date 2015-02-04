#!/bin/bash
#qsub -t 1:50 -q medium_big.q -l mem=6G clscr_exp1_job.sh;
#qsub -t 1:50 -q medium_big.q -l mem=4G clscr_exp2_job.sh;
qsub -q short.q -l mem=6G subgesture_job.sh;