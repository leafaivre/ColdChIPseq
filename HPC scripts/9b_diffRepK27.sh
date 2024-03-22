#!/bin/bash

#SBATCH --job-name=diffreps              # replace name
#SBATCH --mail-user=leafaivre@zedat.fu-berlin.de  # replace email address
#SBATCH --mail-type=end
#SBATCH --nodes=1
#SBATCH --ntasks=12                              # replace with value for your job
#SBATCH --mem-per-cpu=3000                     # replace with value for your job
#SBATCH --time=00:30:00                         # replace with value for your job
#SBATCH --qos=standard                          # replace with value for your job

module load Perl/5.32.0-GCCcore-10.2.0

cd /scratch/leafaivre/ColdChIP_new/results/diffreps               

Con1=$1
Con2=$2
Tre1=$3
Tre2=$4
InCon1=$5
InCon2=$6
InTre1=$7
InTre2=$8
Cond=$9

export PERL5LIB=~/perl5/lib/perl5/:~/perl5/lib/perl5/lib/perl5/site_perl:~/perl5/lib/perl5/lib/perl5/5.32.0

/home/leafaivre/perl5/bin/diffReps.pl \
-tr  "$Tre1".midfrag.bed "$Tre2".midfrag.bed \
-co "$Con1".midfrag.bed "$Con2".midfrag.bed \
--btr "$InTre1".midfrag.bed "$InTre2".midfrag.bed \
--bco "$InCon1".midfrag.bed "$InCon2".midfrag.bed \
-re diff.K27."$Cond".200.midfrag.txt \
-me nb \
--frag 0 \
--nproc 12 \
--noanno \
--nsd broad \
--window 200 \
--chrlen /scratch/leafaivre/ColdChipseq/reference/Chr_Length_TAIR10.txt 


