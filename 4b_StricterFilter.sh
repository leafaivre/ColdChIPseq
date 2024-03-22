#!/bin/bash

#SBATCH --job-name=sorting_mapping              # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1                              # replace with value for your job
#SBATCH --mem-per-cpu=100                     # replace with value for your job
#SBATCH --time=00:30:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load SAMtools/1.9-foss-2018b

cd /scratch/leafaivre/ColdChIP_new/results/filter               

samtools view -b -h -F 4 -q 10 -f 0x2 A${SLURM_ARRAY_TASK_ID}_USR.bam > A${SLURM_ARRAY_TASK_ID}_USR_q10.bam 
samtools index A${SLURM_ARRAY_TASK_ID}_USR_q10.bam

samtools view -c A${SLURM_ARRAY_TASK_ID}_USR_q10.bam > A${SLURM_ARRAY_TASK_ID}.txt
