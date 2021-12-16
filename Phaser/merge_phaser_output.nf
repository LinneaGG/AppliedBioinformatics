imputed_files = Channel.fromPath("/crex/proj/snic2019-8-216/private/albinlinnea/pipelines/phaser/imputed_samples/*postgenos*")

def getLikelihood( file ){
    //regexpPE = /([0-9]\.[0-9])/
    regexpPE = /(\d+\.\d+)/
    (file =~ regexpPE)[0][1]
}

imputed_files
  .map { file -> tuple(getLikelihood(file.baseName), file) }
  .groupTuple()
  .set { grouped_files }

studfile = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/camilles_repository/genotypooler/data/study.population')

imputed_and_studfile = grouped_files.combine(studfile)

process merge_files {
    input:
    tuple prefix, file(imputed_files), stud from imputed_and_studfile
    //file stud from studfile
    //each f from indexChannel

    //worked when having "*vcf.gz*" as output, and !{prefix}.imputed.vcf.gz for merging and sorting
    output:
    file '*merged.imputed.vcf.gz*' into mergedVCF

    publishDir "/crex/proj/snic2019-8-216/private/albinlinnea/pipelines/phaser/merged_imputed_samples/", mode: 'copy'

    shell:
    '''
    #!/bin/bash -l

    for i in !{imputed_files}
    do
            bcftools index -f $i
    done

    #cp $f .
    bcftools merge -Oz -o !{prefix}.merged.imputed.vcf.gz !{imputed_files}

    # sort the merged vcf file
    bcftools view -S !{stud} -Oz -o tmp.sorted.vcf.gz !{prefix}.merged.imputed.vcf.gz
    bcftools view -Oz -o !{prefix}.merged.imputed.vcf.gz tmp.sorted.vcf.gz
    bcftools index -f !{prefix}.merged.imputed.vcf.gz

    '''
}
