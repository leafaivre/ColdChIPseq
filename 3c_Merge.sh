#!/bin/bash

#SBATCH --job-name=merging                # replace name
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

samtools merge -@ 12 "$sample"_merge_unsorted.bam "$sample"_nonfiltered_sorted.bam "$sample"bis_nonfiltered_sorted.bam 

samtools sort -@ 12 -m 4G "$sample"_merge_unsorted.bam -o "$sample"_merge_sorted.bam

samtools index "$sample"_merge_sorted.bam 
