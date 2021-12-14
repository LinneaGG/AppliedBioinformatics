sample = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/IMP.chr20.snps.gt.vcf.gz')
ref =Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/REF.chr20.snps.gt.vcf.gz')
map = Channel.fromPath("/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/plink.chr20.GRCh37.map")
phase = Channel.fromPath("/crex/proj/snic2019-8-216/private/phaser-example/phase")
template_vcf = Channel.fromPath("/crex/proj/snic2019-8-216/private/phaser-example/create_template_vcf.sh")
template_vcf_gtgp = Channel.fromPath("/crex/proj/snic2019-8-216/private/phaser-example/create_template_vcf_gtgp.sh")

csv_script = Channel.fromPath("/home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/likelihood_csv.sh")

process create_csv {
    input:
    file f from csv_script

    output:
    file "likelihoods.csv" into csvChannel

    """
    bash $f -0.75 0.05 -0.4 > likelihoods.csv
    """

}

process likelihood_erros {
    input:
    each likelihood from csvChannel.splitCsv(header:false).flatten()
    file f from sample

    output:
    file "*.vcf" into vcf_channel

    "bash /home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/zipped_blurring_script.sh $likelihood $f > ${f.baseName}_${likelihood}.vcf"

}

process divide_samples {
    time '1h'
    executor 'slurm'
    clusterOptions '-A snic2021-22-462 -p core -n 1 -J divide_samples'

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


process phaser{
    time '10h'
    executor 'slurm'
    clusterOptions '-A snic2021-22-462 -C mem128GB -p node -n 20 -J run_phaser'

    input:
    each s from per_sample_channel.flatten()
    file ref from ref
    file m from map
    file p from phase
    file temp_vcf from template_vcf
    file temp_vcf_gtgp from template_vcf_gtgp

    shell:
    '''
    ne=11418
    error=0.001
    results_directory=/proj/g2021012/nobackup/work/albinlinnea/phaser/imputed_samples/

    # parameters for the linear_state  branch
    algo=integrated
    algo=separate
    algo=fixed
    niter=3

    ./!{temp_vcf} !{s.getParent()}/ !{s.getName()}
    ./!{temp_vcf_gtgp} !{s.getParent()}/ !{s.getName()}
    ./!{p}  --results_directory $results_directory --directory !{s.getParent()}/ --sample_file !{s.getName()} --reference_file !{ref} --ne $ne --map !{m} --error $error

    '''
}
