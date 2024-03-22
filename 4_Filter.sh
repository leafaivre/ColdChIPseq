#!/bin/bash

#SBATCH --job-name=sorting_mapping              # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1                              # replace with value for your job
#SBATCH --mem-per-cpu=8000                     # replace with value for your job
#SBATCH --time=02:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load SAMtools/1.9-foss-2018b

cd /scratch/leafaivre/ColdChIP_new/results               

Sample=$1

samtools view -b -h -F 4 -q 1 align/A${SLURM_ARRAY_TASK_ID}_nonfiltered_sorted.bam > filter/A${SLURM_ARRAY_TASK_ID}_US.bam 
samtools index filter/A${SLURM_ARRAY_TASK_ID}_US.bam

samtools rmdup filter/A${SLURM_ARRAY_TASK_ID}_US.bam filter/A${SLURM_ARRAY_TASK_ID}_USR.bam
samtools index filter/A${SLURM_ARRAY_TASK_ID}_USR.bam

rm align/A${SLURM_ARRAY_TASK_ID}_nonfiltered_sorted.*
rm filter/A${SLURM_ARRAY_TASK_ID}_US.*

