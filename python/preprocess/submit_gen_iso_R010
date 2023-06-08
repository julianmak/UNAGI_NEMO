#!/bin/bash

# NOTE: Lines starting with "#SBATCH" are valid SLURM commands or statements,
#       while those starting with "#" and "##SBATCH" are comments.  Uncomment
#       "##SBATCH" line means to remove one # and start with #SBATCH to be a
#       SLURM command or statement.

#===============================================================
# DEFINE SOME JUNK FOR THE SUBMISSION (??? make this more flexible with e.g. queues?)
#===============================================================

#SBATCH -J iso10            # job name 
#SBATCH -o stdouterr_%j     # output and error file name
#SBATCH -n 40                # total number of mpi tasks requested
#SBATCH -N 1                # total number of nodes requested
#SBATCH -p cpu              # queue (partition) -- standard, development, etc.
##SBATCH -A oces
#SBATCH -t 72:00:00         # maximum runtime
#SBATCH -x hhnode-ib-[201-228] # avoid some nodes

##SBATCH --mail-type=begin        # send email when job begins
##SBATCH --mail-type=end          # send email when job ends
##SBATCH --mail-type=fail         # send email if job fails
##SBATCH --mail-user=jclmak@ust.hk

# Setup runtime environment if necessary
# For example, setup MPI environment
module load openmpi3
python --version

which mpirun

# just make sure the environment really is off otherwise it seems to fail to load the environment here
source /opt/ohpc/pub/apps/anaconda3/bin/activate
source /opt/ohpc/pub/apps/anaconda3/bin/activate py38
python --version
# or you can source ~/.bashrc or ~/.bash_profile

which mpirun

#===============================================================
# LAUNCH JOB
#===============================================================

echo " _ __   ___ _ __ ___   ___         "
echo "| '_ \ / _ \ '_ ' _ \ / _ \        "
echo "| | | |  __/ | | | | | (_) |       "
echo "|_| |_|\___|_| |_| |_|\___/  v3.7  "

echo "OK: ...and here is Christopher doing some calculations for you......"
echo "                  ,-.____,-.          "
echo "                  /   ..   \          "
echo "                 /_        _\         "
echo "                |'o'      'o'|        "
echo "               / ____________ \       "
echo "             , ,'    '--'    '. .     "
echo "            _| |              | |_    "
echo "          /  ' '              ' '  \  "
echo "         (    ',',__________.','    ) "
echo "          \_    ' ._______, '     _/  "
echo "             |                  |     "
echo "             |    ,-.    ,-.    |     "
echo "              \      ).,(      /      "
echo "         gpyy   \___/    \___/        "

# Go to the job submission directory and run your application

# EXP_R010
#export data_dir=/scratch/PI/jclmak/data/users/julian/NEMO/UNAGI/nemo4.0.5/EXP_R010/no_split/gm0000/tau400x/ANALYSIS/
#export data_dir=/scratch/PI/jclmak/data/users/julian/NEMO/UNAGI/nemo4.0.5/EXP_R010/no_split/alp0050/tau100x/ANALYSIS/
#export data_dir=/scratch/PI/jclmak/data/users/julian/NEMO/UNAGI/nemo4.0.5/EXP_R010/split/alp0025/tau100x/ANALYSIS/
#export data_dir=/scratch/PI/jclmak/data/users/julian/NEMO/UNAGI/nemo4.0.5/EXP_R010/split_100km/alp0050/tau100x/ANALYSIS/
#/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/bin/mpirun -n 1 python gen_isodeps_tave.py ${data_dir} "*MOC_T*"

export data_dir=/scratch/PI/jclmak/data/users/julian/NEMO/UNAGI/nemo4.0.5/EXP_R010/split_200km/alp0060_lam80/tau100x/ANALYSIS/
/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.4/bin/mpirun -n 1 python gen_isodeps_tave.py ${data_dir} "*MOC_T*"

wait

sbatch submit_gen_energy_R010

echo "Done!"
