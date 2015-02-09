#!bin/bash
#qsub -q long_big subgesture_job.sh
qsub -t 1 -q long_big.q -l mem=8G -N subgesture subgesture_job.sh;
#for i in {1..10}
#do
#	qsub -t $i -q medium_big.q -N subgesture -hold_jid subgesture subgesture_job.sh;
#done
