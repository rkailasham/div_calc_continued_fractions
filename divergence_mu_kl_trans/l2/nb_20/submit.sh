#!/bin/bash
#SBATCH --job-name=n20
#SBATCH --partition=comp,short
#SBATCH --constraint=Xeon-E5-2667-v3
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=0-10:55:00
#SBATCH --error=n20.out

#module load matlab/r2019b

export MATLABROOT=/usr/local/matlab/r2019b
#mcc -mv driver.m

# The job command(s): 
date=`date +"%H:%M:%S on %d %b %Y"`
echo 
echo "========================================="
echo "Timing: Commenced at $date " 
./run_mu_divergence_quantify_error.sh $MATLABROOT
sleep 5
date=`date +"%H:%M:%S on %d %b %Y"`
echo "Timing: Finished at $date " 
echo "========================================="
