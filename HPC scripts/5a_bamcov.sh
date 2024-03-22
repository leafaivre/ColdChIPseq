#!/bin/bash

#SBATCH --job-name=bamcoverage               # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=12                             
#SBATCH --mem-per-cpu=10                      # replace with value for your job
#SBATCH --time=00:10:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module add deepTools/3.3.1-foss-2020a-Python-3.8.2 

cd /scratch/leafaivre/ColdChIP_new/results/

bamCoverage -p 12 -b filter/A${SLURM_ARRAY_TASK_ID}_USR_q10.bam \
 -o bwtracks/A${SLURM_ARRAY_TASK_ID}_q10.bw \
-bs 10 --skipNonCoveredRegions  \
--normalizeUsing RPKM 
