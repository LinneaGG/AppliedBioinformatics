#!/bin/bash -l
#SBATCH -A snic2021-22-462
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -J plot_beagle_together
# #SBATCH --mail-type=ALL
# #SBATCH --mail-user linnea.gauffingood.6719@student.uu.se

module load bioinfo-tools
module load python/3.6.8

imputed_files="../Beagle_run_samples_together/*.vcf.gz"

truegt="IMP.chr20.snps.gt.vcf.gz"
truegl="truegl.vcf.gz"
outdir="beagle_together"

for imp in $imputed_files
do
        printf "results\n$truegt\n$truegl\n$imp\n$outdir\n5\n0.01\n1" > beagle_argsfile.txt
        python -u eval.py @beagle_argsfile.txt
done
