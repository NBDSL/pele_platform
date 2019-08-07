#!/bin/bash
#SBATCH -J PELE_MPI
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=1000

#############################NO CHANGE###########################
module purge
unset PYTHONPATH
unset LD_LIBRARY_PATH
module load impi 
module load Python/2.7.14-foss-2018a
export PYTHONPATH=/sNow/easybuild/centos/7.4.1708/Skylake/software/PyMOL/2.2.0_0/lib/python2.7/site-packages/:$PYTHONPATH
export PYTHONPATH=/work/NBD_Utilities/PELE/PELE_Softwares/pele_platform/:/work/NBD_Utilities/PELE/PELE_Softwares/pele_platform/pele_platform/:$PYTHONPATH



############################CHANGE##########################

python /work/NBD_Utilities/PELE/PELE_Softwares/pele_platform/tests/test_kinase.py

