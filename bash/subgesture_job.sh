#! /bin/bash

#PBS -N Matlab
#PBS -m be
#PBS -k oe

/opt/matlab/bin/matlab -nodesktop <<EOF
cd('/home/cvc/vponce/Subgesture/');
tic;
runGenDTW;
toc
exit
EOF
