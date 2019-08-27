#! /bin/bash
#$ -S /bin/bash
#$ -V
#$ -pe mpi 4
#$ -cwd

####################################################################
# INFO:
# The script works as following:
# 1/ Set up the case folder whatever you want but not in the scratch
# 2/ Define the path of you scratch folder by setting the SCRATCH variable
# 3/ Submit the job, the calculations will be performed on the scratch
#    and the result will be copied back in the command folder
# 4/ Set the number of cpus you need by defining the entry -pe mpi NCPU
####################################################################

####################################################################
#                                                                                                                                                                           #
#                           Please clean your scratch folder once the calculations are finished                                      #
#                                                                                                                                                                           #
####################################################################

source /u/sw/etc/profile
#module avail
module load gcc-glibc
cd /fast-scratch/mpellegrino/dsmc_parallel
mkdir output_files

# time /home/matematica/barbante/non_ideal_fluid/Programmi/ev_pist_3b.exe &> out
# mpirun -np 4 /home/matematica/barbante/non_ideal_fluid/Programmi/ev_pist_3b.exe &> out
mpiexec /home/matematica/mpellegrino/enskog_vlasov/enskog_vlasov_parallel/main &> main.log

#####################################################################
