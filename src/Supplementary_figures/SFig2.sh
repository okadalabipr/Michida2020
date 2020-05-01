bam="path to your ChIP-seq bam file you are converting into TDF format"
tdf="path to output of this command"

igvtools count -e 200 -w 25 ${bam} ${tdf} mm10