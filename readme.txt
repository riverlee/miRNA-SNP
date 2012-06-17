1) dbsnp135 vcf download at ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/v4.0/00-All.vcf.gz

something to make sure.

1) Is snp position in vcf file 0-based or 1-based

It is 1-based

2) Is pre-miRNA 1-based 
yes


The code is to get SNPs located in miRNA/pre-miRNA and their relationship(in loop,
upstream of miRNA, downstream of miRNA, seed region or just in the miRNA)

miRNA data is from mirBase while dbSNP data is from vcf file


