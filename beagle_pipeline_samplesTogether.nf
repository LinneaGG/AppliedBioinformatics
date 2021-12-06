//nextflow pipeline.nf -c pipeline.config --project "uppmax2020-2-5" --clusterOptions "-M snowy" 

csv_script = Channel.fromPath("/home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/likelihood_csv.sh")

process create_csv {
    input:
    file f from csv_script

    output:
    file "likelihoods.csv" into csvChannel

    publishDir "/home/allu5328/Documents/applied_bioinformatics/Beagle_run_samples_together"

    """
    bash $f -2.0 0.1 -0.3 > likelihoods.csv
    """

}

//likelihood_csv=file("/home/allu5328/Documents/applied_bioinformatics/Beagle_run_samples_together/likelihoods.csv")
//likelihoods = Channel.of(likelihood_csv).splitCsv(header:false).flatten()
study_population = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/IMP.chr20.snps.gt.vcf.gz')
//study_population = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/tests/blurring_script/IMP_subset.vcf')

process likelihood_erros {
    input:
    each likelihood from csvChannel.splitCsv(header:false).flatten()
    file f from study_population

    output:
    file "*.vcf" into vcf_channel

    publishDir "/home/allu5328/Documents/applied_bioinformatics/Beagle_run_samples_together/temp"

    "bash /home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/zipped_blurring_script.sh $likelihood $f > blurred_${likelihood}.vcf"

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

beagle=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/bin/beagle.11Mar19.69c.jar')
reference=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/REF.chr20.snps.gt.vcf.gz')
//reference=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/tests/beagle_test/reference_for_some_markers.vcf')
cfgtjar=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/bin/conform-gt.jar')
beagle_script=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/Beagle_run_samples_together/pipeline_beagle_script.sh')

process beagle {
    time '2h'
    executor 'slurm'
    clusterOptions '-A snic2021-22-462 -C mem128GB -p node -n 1 -J beagle_test --mail-type=ALL --mail-user albin.lundin.5328@student.uu.se'

    input:
    each f from compressed_vcf_channel
    file b from beagle
    file ref from reference
    file cfgt from cfgtjar
    file script from beagle_script

    output:
    file "*blurred.imputed.vcf.gz" into beagle_out
    //file "*imputed*" into outChannel

    publishDir "/home/allu5328/Documents/applied_bioinformatics/Beagle_run_samples_together/real_imputed_vcf"

    /*
    """
    bash $script $b $cfgt $ref $f
    """
    */
    
    shell:
    '''
    #!/bin/bash -l
    
    ## This part was the original code
    ##java -Xss5m -Xmx71680m -jar $b gl=$f \
    ##ref=$ref gprobs=true out=${f.baseName}.imputed

    #beaglejar=~/1000Genomes/scripts/beagle.11Mar19.69c.jar
    #cfgtjar=~/1000Genomes/scripts/conform-gt.jar
    #ref=REF.chr20.snps.gt.vcf.gz
    #vcf=blurred_-1.2.vcf.gz

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




