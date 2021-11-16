csv_script = Channel.fromPath("/home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/likelihood_csv.sh")

process create_csv {
    input:
    file f from csv_script

    output:
    file "likelihoods.csv" into csvChannel

    publishDir "/home/allu5328/Documents/applied_bioinformatics/Nextflow"

    """
    bash $f -2 0.1 0 > likelihoods.csv
    """

}

likelihood_csv=file("/home/allu5328/Documents/applied_bioinformatics/Nextflow/likelihoods.csv")
likelihoods = Channel.of(likelihood_csv).splitCsv(header:false).flatten()
study_population = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/tests/blurring_script/IMP_subset.vcf')

process likelihood_erros {
    input:
    each likelihood from likelihoods
    file f from study_population

    output:
    file "*.vcf" into vcf_channel

    "bash /home/allu5328/Documents/applied_bioinformatics/AppliedBioinformatics/blurring_script.sh $likelihood $f > test_${likelihood}.vcf"

}

beagle=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/bin/beagle.11Mar19.69c.jar')
reference=Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/tests/beagle_test/reference_for_some_markers.vcf')

process beagle {
    input:
    each f from vcf_channel
    file b from beagle
    file ref from reference

    output:
    file "*" into outChannel

    publishDir "/home/allu5328/Documents/applied_bioinformatics/Nextflow/vcf_files"

    """
    java -Xss5m -jar $b gl=$f \
    ref=$ref impute=true gprobs=true out=${f.baseName}.imputed
    """
}

outChannel.view()
