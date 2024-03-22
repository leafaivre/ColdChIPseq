#!/bin/bash

#SBATCH --job-name=Matrix                   # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=12                             # replace with value for your job
#SBATCH --mem-per-cpu=30                      # replace with value for your job
#SBATCH --time=00:10:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load deepTools

cd /scratch/leafaivre/ColdChIP_new/results

List=$1
Modif=$2
Name=$3

computeMatrix scale-regions \
-b 500 -a 500 \
--regionBodyLength 2000 \
-R bedFiles/"$List".bed \
-S bwtracks/N_"$Modif"_merged.bw \
bwtracks/3h_"$Modif"_merged.bw \
bwtracks/3d_"$Modif"_merged.bw \
--skipZeros \
-o metagene/"$Name".gz \
-p 12


plotProfile -m metagene/"$Name".gz -out metagene/"$Name".png --perGroup \
 --plotTitle "" --samplesLabel "N" "3h" "3d" \
 --colors green skyblue darkblue \
 --refPointLabel "TSS" -T "" -z ""

