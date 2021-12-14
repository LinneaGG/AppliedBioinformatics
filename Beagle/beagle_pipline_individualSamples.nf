csv_script = Channel.fromPath("/home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/likelihood_csv.sh")

process create_csv {
    input:
    file f from csv_script

    output:
    file "likelihoods.csv" into csvChannel

    """
    bash $f -0.8 0.1 -0.4 > likelihoods.csv
    """

}

study_population = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/IMP.chr20.snps.gt.vcf.gz')

process likelihood_erros {
    input:
    each likelihood from csvChannel.splitCsv(header:false).flatten()
    file f from study_population

    output:
    file "*.vcf" into vcf_channel

    "bash /home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/zipped_blurring_script.sh $likelihood $f > ${f.baseName}_${likelihood}.vcf"

}

process divide_samples {
    time '1h'
    executor 'slurm'
    clusterOptions '-A snic2021-22-462 -p core -n 1 -J beagle_run_per_sample'

    input:
    file f from vcf_channel

    output:
    file "*.vcf.gz" into per_sample_channel

    shell:
    '''
    #!/bin/bash -l

    tot_cols=$(awk '{print NF}' !{f} | sort -nu | tail -n 1)
    sample_cols=$(expr $tot_cols - 9)

    for i in $(seq 1 1 $sample_cols)
    do
        sampleId=$( bcftools query -l !{f} | head -$i | tail -1 )
        sample_file=s$sampleId.!{f}.gz
        bcftools view -Oz -s $sampleId !{f} > $sample_file
    done
    '''
}



beagle=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/bin/beagle.11Mar19.69c.jar')
reference=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/REF.chr20.snps.gt.vcf.gz')
cfgtjar=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/bin/conform-gt.jar')

process beagle {
    time '1h'
    executor 'slurm'
    clusterOptions '-A snic2021-22-462 -p core -n 1 -J beagle_run_per_sample'


    input:
    each f from per_sample_channel.flatten()
    file b from beagle
    file ref from reference
    file cfgt from cfgtjar

    output:
    file "*blurred.imputed.vcf.gz" into beagleOutChannel

    shell:
    '''

    chrom=$( bcftools query -f '%CHROM\n' !{ref} | head -1)
    startpos=$( bcftools query -f '%POS\n' !{ref} | head -1 )
    endpos=$( bcftools query -f '%POS\n' !{ref} | tail -1 )

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


def getLikelihood( file ){
    //regexpPE = /([0-9]\.[0-9])/
    regexpPE = /(\d+\.\d+)/
    (file =~ regexpPE)[0][1]
}

beagleOutChannel
  .map { file -> tuple(getLikelihood(file.baseName), file) }
  .groupTuple()
  .set { grouped_files }

studfile = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/study.population')

beagle_and_studfile = grouped_files.combine(studfile)

process merge_files {
    input:
    tuple prefix, beagleFiles, stud from beagle_and_studfile

    output:
    file '*merged*imputed.vcf.gz*' into mergedVCF

    publishDir "/crex/proj/snic2019-8-216/private/albinlinnea/pipelines/beagle_individually/merged_vcfs", mode: 'copy'

    shell:
    '''
    #!/bin/bash -l

    for i in !{beagleFiles}
    do
            bcftools index $i
    done

    #cp $f .
    bcftools merge -Oz -o merged.-!{prefix}.imputed.vcf.gz !{beagleFiles.join(' ')}

    # sort the merged vcf file
    bcftools view -S !{stud} -Oz -o tmp.sorted.vcf.gz merged.-!{prefix}.imputed.vcf.gz
    bcftools view -Oz -o merged.-!{prefix}.imputed.vcf.gz tmp.sorted.vcf.gz
    bcftools index -f merged.-!{prefix}.imputed.vcf.gz

    '''
}
