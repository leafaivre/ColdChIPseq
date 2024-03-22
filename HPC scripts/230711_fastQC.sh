#!/bin/bash

#SBATCH --job-name=fastQC_A1to24               # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=500                      # replace with value for your job
#SBATCH --time=00:30:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module add FastQC/0.11.7-Java-1.8.0_162
         # replace with value for your job


fastqc -o /scratch/leafaivre/ColdChIP_new/results/fastqc /scratch/leafaivre/ColdChipseq/raw_data/A${SLURM_ARRAY_TASK_ID}/*.fq.gz
