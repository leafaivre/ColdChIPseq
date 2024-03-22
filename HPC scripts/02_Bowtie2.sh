#!/bin/bash

#SBATCH --job-name=alignment                   # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=12                             # replace with value for your job
#SBATCH --mem-per-cpu=200                      # replace with value for your job
#SBATCH --time=05:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module add Bowtie2

cd /scratch/leafaivre

sample=$1

bowtie2 -p 6 --local --no-unal  -x ColdChipseq/reference/TAIR10/TAIR10 \
-1 ColdChipseq/raw_data/"$sample"/"$sample"_*_1.fq.gz -2 ColdChipseq/raw_data/"$sample"/"$sample"_*_2.fq.gz \
-S ColdChIP_new/results/align/"$sample"_nonfiltered_NoUnal.sam \
2> ColdChIP_new/results/align/"$sample".log

