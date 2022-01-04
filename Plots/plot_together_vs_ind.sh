#!/bin/bash -l
#SBATCH -A uppmax2020-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -J plot_beagle_together
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user linnea.gauffingood.6719@student.uu.se

module load bioinfo-tools
module load python/3.6.8

imp1=(../Results/BeagleTogether/*.vcf.gz)
imp2=(../Results/BeagleIndividually/*.vcf.gz)

truegt="../data/IMP.chr20.snps.gt.vcf.gz"
truegl="truegl.vcf.gz"
outdir="beagle_together_vs_ind"

for ((i=0;i<${#imp1[@]};i++)); do
   printf "results\n$truegt\n$truegl\n${imp1[i]}\n${imp2[i]}\n$outdir\n5\n0.01\n1" > two_argsfile.txt
        python -u eval_together_vs_ind.py @two_argsfile.txt

done
