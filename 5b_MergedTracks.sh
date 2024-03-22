#!/bin/bash

#SBATCH --job-name=merging                # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=12                             
#SBATCH --mem-per-cpu=4000                      # replace with value for your job
#SBATCH --time=02:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load SAMtools/1.9-foss-2018b

cd /scratch/leafaivre/ColdChIP_new/results/

sample1=$1
sample2=$2
Cond=$3

samtools merge -@ 12  bwtracks/"$Cond"_unsorted.bam filter/"$sample1"_USR_q10.bam filter/"$sample2"_USR_q10.bam

samtools sort -@ 12 -m 4G bwtracks/"$Cond"_unsorted.bam -o bwtracks/"$Cond"_sorted.bam

samtools index bwtracks/"$Cond"_sorted.bam 

rm bwtracks/"$Cond"_unsorted.bam

module add deepTools/3.3.1-foss-2020a-Python-3.8.2 

cd bwtracks/

bamCoverage -p 12 -b "$Cond"_sorted.bam \
 -o "$Cond"_merged.bw \
-bs 10 --skipNonCoveredRegions \
--normalizeUsing RPKM 
