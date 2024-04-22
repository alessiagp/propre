#!/bin/bash

#################################################################################################################################
### ATTENTION: Hereafter, an example for executing MPI scripts on 48 cores on UniTn Cluster. Change it for your cluster/machine 

#PBS -l select=1:ncpus=24:mpiprocs=1:mem=4gb
#PBS -q VARIAMOLS_cpuQ

echo "loading modules"
module load mpich-3.2

echo "modules loaded"

echo "running simulations"

cd $PBS_O_WORKDIR
#################################################################################################################################

# Executing serial version of "block.py" and then "CANVAS.py" for the APO-2 conformation of Pembrolizumab Antibody 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-2.txt, that will be the input of CANVAS.py script  

name="apo3-2A_noH"

PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/${name}


python3 $PYTHONDIR/ResRel-MPI.py -r $inputDIR/${name}.gro -t $inputDIR/${name}.xtc  

rm -r test-${name}
mkdir test-${name}

mv trace_${name}.txt Hs-Hk-Nsites-${name}.txt test-${name}
