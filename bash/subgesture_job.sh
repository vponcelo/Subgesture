#! /bin/bash

#PBS -N Matlab
#PBS -m be
#PBS -k oe

/opt/matlab/bin/matlab -nodesktop <<EOF
cd('/home/cvc/vponce/Subgesture/');
tic
runGenDTW('overlap',$SGE_TASK_ID);
toc
exit
EOF
