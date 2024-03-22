#!/bin/bash

#SBATCH --job-name=IDR                # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=1                             
#SBATCH --mem-per-cpu=500                      # replace with value for your job
#SBATCH --time=00:35:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load Anaconda3/2020.02
cd ~/idr-2.0.3
python3 setup.py install --user


cd /scratch/leafaivre/ColdChIP_new/results/idr

sample1=$1
sample2=$2

python3 ~/idr-2.0.2/bin/idr --samples "$sample1"_sortedpeaks.broadPeak "$sample2"_sortedpeaks.broadPeak \
--input-file-type broadPeak \
--rank p.value \
--output-file "$sample1""$sample2"-idr \
--plot \
--log-output-file "$sample1""$sample2"-idr.log

