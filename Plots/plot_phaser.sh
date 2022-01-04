#!/bin/bash -l
#SBATCH -A uppmax2020-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J plot_phaser
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user linnea.gauffingood.6719@student.uu.se

module load bioinfo-tools
module load python/3.6.8

imputed_files="../Results/Phaser_merged/*.vcf.gz"

truegt="../data/IMP.chr20.snps.gt.vcf.gz"
truegl="truegl.vcf.gz"
outdir="phaser_220103"

for imp in $imputed_files
do
        printf "results\n$truegt\n$truegl\n$imp\n$outdir\n5\n0.01\n1" > phaser_argsfile.txt
	python -u eval_phaser.py @phaser_argsfile.txt
done


