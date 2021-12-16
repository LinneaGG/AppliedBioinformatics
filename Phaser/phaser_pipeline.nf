filePath="$PWD"

sample = Channel.fromPath(filePath + '/../data/IMP.chr20.snps.gt.vcf.gz')
ref =Channel.fromPath(filePath + '/../data/REF.chr20.snps.gt.vcf.gz')
map = Channel.fromPath(filePath + "/../data/plink.chr20.GRCh37.map")
phase = Channel.fromPath(filePath + "/../bin/Phaser/phase")
template_vcf = Channel.fromPath(filePath + "/../bin/Phaser/create_template_vcf.sh")
template_vcf_gtgp = Channel.fromPath(filePath + "/../bin/Phaser/create_template_vcf_gtgp.sh")

csv_script = Channel.fromPath(filePath + "/../Scripts/likelihood_csv.sh")

process create_csv {
    input:
    file f from csv_script

    output:
    file "likelihoods.csv" into csvChannel

    """
    bash $f -0.75 0.05 -0.75 > likelihoods.csv
    """

}

process likelihood_erros {
    input:
    each likelihood from csvChannel.splitCsv(header:false).flatten()
    file f from sample

    output:
    file "*.vcf" into vcf_channel

    "bash ${filePath}/../Scripts/zipped_blurring_script.sh $likelihood $f > ${f.baseName}_${likelihood}.vcf"

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
    results_directory=!{filePath}/../Results/Phaser_imputed/

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
