//nextflow pipeline.nf -c pipeline.config --project "uppmax2020-2-5" --clusterOptions "-M snowy" 

csv_script = Channel.fromPath("/home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/likelihood_csv.sh")

process create_csv {
    input:
    file f from csv_script

    output:
    file "likelihoods.csv" into csvChannel

    //publishDir "/home/allu5328/Documents/applied_bioinformatics/Nextflow"

    """
    bash $f -2.0 0.1 -0.4 > likelihoods.csv
    """

}

//likelihood_csv=file("/home/allu5328/Documents/applied_bioinformatics/Nextflow/likelihoods.csv")
//likelihoods = Channel.of(likelihood_csv).splitCsv(header:false).flatten()
study_population = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/IMP.chr20.snps.gt.vcf.gz')
//study_population = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/Nextflow/IMP_subset_lessSamples.vcf')
//study_population = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/tests/blurring_script/IMP_subset.vcf')

process likelihood_erros {
    input:
    each likelihood from csvChannel.splitCsv(header:false).flatten()
    file f from study_population

    output:
    file "*.vcf" into vcf_channel

    //publishDir "/home/allu5328/Documents/applied_bioinformatics/Nextflow/blurred_vcf"

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

    publishDir "/home/allu5328/Documents/applied_bioinformatics/Nextflow/per_sample_vcf_files"

    shell:
    '''
    #!/bin/bash -l

    #samples_file=reference_for_some_markers.vcf
    #results_directory=samples/
    
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

// probably a better way of creating a file for all samples
//sampleId=$( bcftools query -l $indir$samples_file | head -$sample | tail -1 )
//sample_file=s$sampleId.$samples_file
//bcftools view -Oz -s $sampleId -o $indir$sample_file $indir$samples_file



beagle=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/bin/beagle.11Mar19.69c.jar')
reference=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/REF.chr20.snps.gt.vcf.gz')
//reference=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/tests/beagle_test/reference_for_some_markers.vcf')
//reference=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/tests/beagle_test/reference_for_some_markers.vcf')
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
    //file "*csi" into indexChannel

    //publishDir "/home/allu5328/Documents/applied_bioinformatics/Nextflow/vcf_files"

    /*
    # this is original code that worked for test data{
    #java -Xss5m -jar $b gl=$f \
    #ref=$ref impute=true gprobs=true out=${f.baseName}.imputed
    #bcftools index ${f.baseName}.imputed.vcf.gz (does not need)
    #}

    */

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

process merge_files {
    input:
    set prefix, file(beagleOutChannel) from grouped_files

    output:
    file '*vcf.gz' into mergedVCF

    publishDir "/home/allu5328/Documents/applied_bioinformatics/Beagle_run_per_sample/merged_vcfs"

    shell:
    '''
    #!/bin/bash -l

    for i in !{beagleOutChannel}
    do
	    bcftools index $i
    done

    bcftools merge -Oz -o !{prefix}.imputed.vcf.gz !{beagleOutChannel}
    '''
    //merge_file.sh ${edited_files} -o $prefix.FINAL.bf
}
