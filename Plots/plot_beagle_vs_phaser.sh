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

# Folders must contain the same uncertainties 

imp_beagle=(../Results/BeagleIndividually/*.vcf.gz)
imp_phaser=(../Results/Phaser_merged/*.vcf.gz)

truegt="IMP.chr20.snps.gt.vcf.gz"
truegl="truegl.vcf.gz"
outdir="beagle_vs_phaser"

for ((i=0;i<${#imp_beagle[@]};i++)); do
   	printf "results\n$truegt\n$truegl\n${imp_beagle[i]}\n${imp_phaser[i]}\n$outdir\n5\n0.01\n1" > bp_argsfile.txt
       	python -u eval_beagle_vs_phaser.py @bp_argsfile.txt

done
