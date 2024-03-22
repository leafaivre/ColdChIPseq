#!/bin/bash

#SBATCH --job-name=MACS2               # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --ntasks=1                            
#SBATCH --mem-per-cpu=2000                      # replace with value for your job
#SBATCH --time=00:40:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load MACS2/2.2.5-foss-2018b-Python-3.6.6

cd /scratch/leafaivre/ColdChIP_new/results/

sample=$1

input=$2

macs2 callpeak -t filter/"$sample"_USR_q10.bam \
		-c filter/"$input"_USR_q10.bam \
		-f BAMPE \
		-g 119481543 \
		-n "$sample"_USR \
		-B \
		--broad \
		--keep-dup 1 \
		-p 0.01 \
		--outdir macs2/ \
		2> macs2/"$sample"USR.log





