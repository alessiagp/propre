#!/bin/bash

#################################################################################################################################
### ATTENTION: Hereafter, an example for executing MPI scripts on 48 cores on UniTn Cluster. Change it for your cluster/machine 

#PBS -l select=1:ncpus=36:mpiprocs=1:mem=4gb
#PBS -q VARIAMOLS_cpuQ

echo "loading modules"
module load mpich-3.2

echo "modules loaded"

echo "running simulations"

cd $PBS_O_WORKDIR
#################################################################################################################################

# Executing serial version of "block.py" and then "CANVAS.py" for the APO-2 conformation of Pembrolizumab Antibody 
# NOTE: The output of block.py is the list of survived atoms, i.e. list-atoms-opt-2.txt, that will be the input of CANVAS.py script  

name="ake"

PYTHONDIR=../PYTHON-scripts

inputDIR=../input-files/${name}


## 1st part: remove_H_atoms.py

python3 $PYTHONDIR/remove_H_atoms.py -r $inputDIR/${name}.gro -t $inputDIR/${name}.xtc

mv Reference_noH.gro ../input-files/${name}/${name}_noH.gro 
mv Trajectory_noH.xtc ../input-files/${name}/${name}_noH.xtc

#############################################################

## 2nd part: ResRel.py

name_noH="${name}_noH"

python3 $PYTHONDIR/ResRel-MPI.py -r $inputDIR/${name_noH}.gro -t $inputDIR/${name_noH}.xtc  

rm -r test-${name_noH}
mkdir test-${name_noH}

mv trace_${name_noH}.txt Hs-Hk-Nsites-${name_noH}.txt test-${name_noH}
