#!/bin/bash
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de
#SBATCH --job-name=ChIQC
#SBATCH --ntasks=12                             
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4096
#SBATCH --time=04:00:00
#SBATCH --qos=standard

module load R/4.0.0-foss-2020a
module load librsvg/2.50.2-foss-2020a

cd /scratch/leafaivre/ColdChIP_new/results

Rscript ~/ColdChIP_new/7_ChIPQC.R 
