#!/bin/bash

#SBATCH --job-name=SAMtoBAM                # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=12                             
#SBATCH --mem-per-cpu=5000                      # replace with value for your job
#SBATCH --time=02:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load SAMtools/1.9-foss-2018b

cd /scratch/leafaivre/ColdChIP_new/results/align

sample=$1

samtools view -h -S -b -o A${SLURM_ARRAY_TASK_ID}bis_nonfiltered_unsorted.bam A${SLURM_ARRAY_TASK_ID}bis_nonfiltered_NoUnal.sam

samtools sort -@ 12 -m 4G A${SLURM_ARRAY_TASK_ID}bis_nonfiltered_unsorted.bam -o A${SLURM_ARRAY_TASK_ID}bis_nonfiltered_sorted.bam

samtools index A${SLURM_ARRAY_TASK_ID}bis_nonfiltered_sorted.bam 


rm A${SLURM_ARRAY_TASK_ID}bis_nonfiltered_NoUnal.sam
rm A${SLURM_ARRAY_TASK_ID}bis_nonfiltered_unsorted.bam
