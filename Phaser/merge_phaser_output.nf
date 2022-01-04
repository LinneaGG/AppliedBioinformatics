filePath="$PWD"

imputed_files = Channel.fromPath(filePath + "/../Results/Phaser_imputed/*postgenos*")

def getLikelihood( file ){
    //regexpPE = /([0-9]\.[0-9])/
    regexpPE = /(\d+\.\d+)/
    (file =~ regexpPE)[0][1]
}

imputed_files
  .map { file -> tuple(getLikelihood(file.baseName), file) }
  .groupTuple()
  .set { grouped_files }

studfile = Channel.fromPath(filePath + '/../data/study.population')

imputed_and_studfile = grouped_files.combine(studfile)

process merge_files {
    input:
    tuple prefix, file(imputed_files), stud from imputed_and_studfile
    //file stud from studfile
    //each f from indexChannel

    //worked when having "*vcf.gz*" as output, and !{prefix}.imputed.vcf.gz for merging and sorting
    output:
    file '*merged.imputed.vcf.gz*' into mergedVCF

    publishDir "${filePath}/../Results/Phaser_merged/", mode: 'copy'

    shell:
    '''
    #!/bin/bash -l

    for i in !{imputed_files}
    do
            bcftools index -f $i
    done

    #cp $f .
    bcftools merge -Oz -o  phaser_-!{prefix}.merged.imputed.vcf.gz !{imputed_files}

    # sort the merged vcf file
    bcftools view -S !{stud} -Oz -o tmp.sorted.vcf.gz  phaser_-!{prefix}.merged.imputed.vcf.gz
    bcftools view -Oz -o  phaser_-!{prefix}.merged.imputed.vcf.gz tmp.sorted.vcf.gz
    bcftools index -f  phaser_-!{prefix}.merged.imputed.vcf.gz

    '''
}
