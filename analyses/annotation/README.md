# Folder with variant annotation scripts for the project

## 1) commands to prepare the dbNSFP annotation file

### create a folder for the annotation files and go there
```
mkdir -p shared/annotation/hg38/dbnsfp
cd shared/annotation/hg38/dbnsfp
```

### download dbNSFP
```
wget https://dbnsfp.s3.amazonaws.com/dbNSFP4.4a.zip
```

### unzip
```
unzip dbNSFP4.4a.zip
```

### extract the header from the first file
```
zcat dbNSFP4.4a_variant.chr1.gz | head -n1 > header.txt
```

### concatenate all files, remove header and add the header from the first file
```
(cat header.txt; zcat dbNSFP4.4a_variant.chr1.gz dbNSFP4.4a_variant.chr2.gz dbNSFP4.4a_variant.chr3.gz dbNSFP4.4a_variant.chr4.gz dbNSFP4.4a_variant.chr5.gz dbNSFP4.4a_variant.chr6.gz dbNSFP4.4a_variant.chr7.gz dbNSFP4.4a_variant.chr8.gz dbNSFP4.4a_variant.chr9.gz dbNSFP4.4a_variant.chr10.gz dbNSFP4.4a_variant.chr11.gz dbNSFP4.4a_variant.chr12.gz dbNSFP4.4a_variant.chr13.gz dbNSFP4.4a_variant.chr14.gz dbNSFP4.4a_variant.chr15.gz dbNSFP4.4a_variant.chr16.gz dbNSFP4.4a_variant.chr17.gz dbNSFP4.4a_variant.chr18.gz dbNSFP4.4a_variant.chr19.gz dbNSFP4.4a_variant.chr20.gz dbNSFP4.4a_variant.chr21.gz dbNSFP4.4a_variant.chr22.gz dbNSFP4.4a_variant.chrX.gz dbNSFP4.4a_variant.chrY.gz dbNSFP4.4a_variant.chrM.gz | grep -v "^#" ) | bgzip -c > dbNSFP4.4a_grch38.gz
```

### index the file with tabix
```
tabix -s 1 -b 2 -e 2 dbNSFP4.4a_grch38.gz
```

### remove the intermediate files
```
rm dbNSFP4.4a_variant.chr*.gz
rm h dbNSFP4.4a.readme.txt search_dbNSFP44a.jar search_dbNSFP44a.readme.pdf search_dbNSFP44a.class tryhg19.in tryhg38.in LICENSE.txt tryhg18.in try.vcf
```