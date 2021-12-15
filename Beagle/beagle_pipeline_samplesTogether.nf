filePath = "$PWD"

csv_script = Channel.fromPath(filePath + "/../Scripts/likelihood_csv.sh")

process create_csv {
    input:
    file f from csv_script

    output:
    file "likelihoods.csv" into csvChannel

    """
    bash $f -2.0 0.1 -2.0 > likelihoods.csv
    """

}

study_population = Channel.fromPath(filePath + '/../data/IMP.chr20.snps.gt.vcf.gz')

process likelihood_erros {
    input:
    each likelihood from csvChannel.splitCsv(header:false).flatten()
    file f from study_population

    output:
    file "*.vcf" into vcf_channel

    "bash ${filePath}/../Scripts/zipped_blurring_script.sh $likelihood $f > blurred_${likelihood}.vcf"

}

process bgzip {
    input:
    file f from vcf_channel

    output:
    file "*.vcf.gz" into compressed_vcf_channel

    """
    bcftools view $f -Oz -o ${f}.gz
    bcftools index ${f}.gz
    """

}

beagle=Channel.fromPath(filePath + '/../bin/beagle.11Mar19.69c.jar')
reference=Channel.fromPath(filePath + '/../data/REF.chr20.snps.gt.vcf.gz')
cfgtjar=Channel.fromPath(filePath + '/../bin/conform-gt.jar')
//beagle_script=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/Beagle_run_samples_together/pipeline_beagle_script.sh')

process beagle {
    time '2h'
    executor 'slurm'
    clusterOptions '-A snic2021-22-462 -C mem128GB -p node -n 1 -J beagle_test --mail-type=ALL --mail-user albin.lundin.5328@student.uu.se'

    input:
    each f from compressed_vcf_channel
    file b from beagle
    file ref from reference
    file cfgt from cfgtjar
//    file script from beagle_script

    output:
    file "*blurred.imputed.vcf.gz" into beagle_out

    publishDir "${filePath}/../Results/BeagleTogether/"
    
    shell:
    '''
    #!/bin/bash -l

    chrom=$( bcftools query -f '%CHROM\n' !{ref} | head -1)
    startpos=$( bcftools query -f '%POS\n' !{f} | head -1 )
    endpos=$( bcftools query -f '%POS\n' !{f} | tail -1 )

    # phasing round 1
    java -Xss5m -jar !{b} impute=false gtgl=!{f} out=!{f.baseName}_unphased.blurred

    # phasing round 2 + sort
    java -Xss5m -jar !{b} impute=false gt=!{f.baseName}_unphased.blurred.vcf.gz out=!{f.baseName}_phased.blurred
    bcftools index -f !{f.baseName}_phased.blurred.vcf.gz
    bcftools sort -Oz -o !{f.baseName}_sorted.phased.blurred.vcf.gz !{f.baseName}_phased.blurred.vcf.gz
    bcftools index -f !{f.baseName}_sorted.phased.blurred.vcf.gz

    # dedup
    bcftools norm --rm-dup all -Oz -o !{f.baseName}_phased.dedup.blurred.vcf.gz !{f.baseName}_sorted.phased.blurred.vcf.gz
    bcftools index -f !{f.baseName}_phased.dedup.blurred.vcf.gz

    # conform-gt
    java -jar !{cfgt} ref=!{ref} gt=!{f.baseName}_phased.dedup.blurred.vcf.gz chrom=$chrom:$startpos-$endpos out=!{f.baseName}_cfgt.blurred
    
    bcftools index -f !{f.baseName}_cfgt.blurred.vcf.gz
    bcftools index -f !{ref}

    # imputation
    java -Xss5m -jar !{b} impute=true gt=!{f.baseName}_cfgt.blurred.vcf.gz ref=!{ref} gprobs=true out=!{f.baseName}_blurred.imputed
    bcftools index -f !{f.baseName}_blurred.imputed.vcf.gz
    '''
    
}




