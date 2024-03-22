#!/bin/bash

#SBATCH --job-name=sortpeaks                # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=500                      # replace with value for your job
#SBATCH --time=00:35:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job


cd /scratch/leafaivre/ColdChIP_new/results/

sample=${SLURM_ARRAY_TASK_ID}

sort -k8,8nr macs2/A"$sample"_USR_peaks.broadPeak > idr/A"$sample"_sortedpeaks.broadPeak


