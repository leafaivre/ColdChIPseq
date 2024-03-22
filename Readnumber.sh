#!/bin/bash

#SBATCH --job-name=sorting_mapping              # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1                              # replace with value for your job
#SBATCH --mem-per-cpu=100                     # replace with value for your job
#SBATCH --time=00:10:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load SAMtools/1.9-foss-2018b

cd /scratch/leafaivre/ColdChIP_new/results/filter               

samtools view -c A${SLURM_ARRAY_TASK_ID}_USR.bam > A${SLURM_ARRAY_TASK_ID}_all.txt
samtools view -c -q 5 A${SLURM_ARRAY_TASK_ID}_USR.bam > A${SLURM_ARRAY_TASK_ID}_5.txt
samtools view -c -q 10 A${SLURM_ARRAY_TASK_ID}_USR.bam > A${SLURM_ARRAY_TASK_ID}_10.txt
