#!bin/bash
qsub -t 1 -q short.q -l mem=2G -N hmm hmmGesture_job.sh;
