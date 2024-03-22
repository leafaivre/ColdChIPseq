#!/bin/bash

#SBATCH --job-name=BedpeToBed6              # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=12                              # replace with value for your job
#SBATCH --mem-per-cpu=2000                     # replace with value for your job
#SBATCH --time=03:00:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load SAMtools/1.9-foss-2018b
module load BEDTools/2.27.1-foss-2018b

cd /scratch/leafaivre/ColdChIP_new/results               

sample=${SLURM_ARRAY_TASK_ID}

samtools sort -n -@ 12 -m 4G filter/A"$sample"_USR_q10.bam -o diffreps/A"$sample"_USR_q10_nsorted.bam

cd diffreps

bedtools bamtobed -bedpe -i A"$sample"_USR_q10_nsorted.bam > A"$sample"_USR_q10.bedpe     

rm A"$sample"_USR_q10_nsorted.bam    

cat A"$sample"_USR_q10.bedpe |\
awk '{chr=$1; start=$2; end=$6; mid=(start+end)/2; printf "%s\t%.1f\t%.1f\t%s\t%d\t.\n", chr,mid,mid+1,$7,$8}' |\
sed 's/\.0//g;s/\.5//g' > A"$sample".midfrag.bed

rm A"$sample"_USR_q10.bedpe
