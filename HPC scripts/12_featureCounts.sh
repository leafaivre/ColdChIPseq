#!/bin/bash

#SBATCH --job-name=featureCounts               # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --ntasks=12
#SBATCH --nodes=1                             
#SBATCH --mem-per-cpu=500                      # replace with value for your job
#SBATCH --time=01:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load Subread/2.0.1-GCC-9.3.0

cd /scratch/leafaivre/ColdChIP_new/results/

featureCounts -T 12 \
-a /scratch/leafaivre/References/TAIR10_noprom.saf \
-F SAF \
-o featureCounts/FC_all_noprom.txt \
-p \
-Q 10 \
filter/*_USR_q10.bam

featureCounts -T 12 \
-a /scratch/leafaivre/References/TAIR10_prom.saf \
-F SAF \
-o featureCounts/FC_all_prom.txt \
-p \
-Q 10 \
filter/*_USR_q10.bam




