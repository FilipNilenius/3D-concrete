#!/usr/bin/env bash
#SBATCH -A C3SE2017-1-4
#SBATCH -p hebbe
#SBATCH -J JobTest
###SBATCH -C "MEM512|MEM1024"
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -o output.stdout
#SBATCH -e output.stderr


### Line for MATLAB with you input file(s)
#################################################################
module load MATLAB
matlab -nodesktop -nosplash -r main
################################################################
 
# End script