csvfile=file("err.csv")
errorChannel=Channel.of(csvfile).splitCsv(header:false).flatten()
fileChannel = Channel.fromPath('/home/allu5328/Documents/applied_bioinformatics/blurring_script/IMP_subset.vcf')

process test{
input:
each x from errorChannel
file f from fileChannel

output:
file "*.vcf" into outChannel

publishDir "/home/allu5328/Documents/applied_bioinformatics/blurring_script"

"bash /home/allu5328/Documents/applied_bioinformatics/blurring_script/script.sh $x $f > test_${x}.vcf"
}

outChannel.view()


