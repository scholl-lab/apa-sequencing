# Folder with variant annotation scripts for the project

## 1) commands to prepare the dbNSFP annotation file

### create a folder for the annotation files and go there
```
mkdir -p shared/annotation/hg38/dbnsfp
cd shared/annotation/hg38/dbnsfp
```

### download dbNSFP
```
# version dbNSFP4.4a
wget https://dbnsfp.s3.amazonaws.com/dbNSFP4.4a.zip

# version dbNSFP4.7a (2024-05-16)
wget https://dbnsfp.s3.amazonaws.com/dbNSFP4.7a.zip
```

### unzip
```
# version dbNSFP4.4a
unzip dbNSFP4.4a.zip

# version dbNSFP4.7a (2024-05-16)
unzip dbNSFP4.7a.zip
```

### extract the header from the first file
```
# version dbNSFP4.4a
zcat dbNSFP4.4a_variant.chr1.gz | head -n1 > header.txt

# version dbNSFP4.7a (2024-05-16)
zcat dbNSFP4.7a_variant.chr1.gz | head -n1 > dbNSFP4.7a_header.txt
```

### concatenate all files, remove header and add the header from the first file
```
# version dbNSFP4.4a
(cat header.txt; zcat dbNSFP4.4a_variant.chr1.gz dbNSFP4.4a_variant.chr2.gz dbNSFP4.4a_variant.chr3.gz dbNSFP4.4a_variant.chr4.gz dbNSFP4.4a_variant.chr5.gz dbNSFP4.4a_variant.chr6.gz dbNSFP4.4a_variant.chr7.gz dbNSFP4.4a_variant.chr8.gz dbNSFP4.4a_variant.chr9.gz dbNSFP4.4a_variant.chr10.gz dbNSFP4.4a_variant.chr11.gz dbNSFP4.4a_variant.chr12.gz dbNSFP4.4a_variant.chr13.gz dbNSFP4.4a_variant.chr14.gz dbNSFP4.4a_variant.chr15.gz dbNSFP4.4a_variant.chr16.gz dbNSFP4.4a_variant.chr17.gz dbNSFP4.4a_variant.chr18.gz dbNSFP4.4a_variant.chr19.gz dbNSFP4.4a_variant.chr20.gz dbNSFP4.4a_variant.chr21.gz dbNSFP4.4a_variant.chr22.gz dbNSFP4.4a_variant.chrX.gz dbNSFP4.4a_variant.chrY.gz dbNSFP4.4a_variant.chrM.gz | grep -v "^#" ) | bgzip -c > dbNSFP4.4a_grch38.gz

# version dbNSFP4.7a (2024-05-16)
(cat dbNSFP4.7a_header.txt; zcat dbNSFP4.7a_variant.chr1.gz dbNSFP4.7a_variant.chr2.gz dbNSFP4.7a_variant.chr3.gz dbNSFP4.7a_variant.chr4.gz dbNSFP4.7a_variant.chr5.gz dbNSFP4.7a_variant.chr6.gz dbNSFP4.7a_variant.chr7.gz dbNSFP4.7a_variant.chr8.gz dbNSFP4.7a_variant.chr9.gz dbNSFP4.7a_variant.chr10.gz dbNSFP4.7a_variant.chr11.gz dbNSFP4.7a_variant.chr12.gz dbNSFP4.7a_variant.chr13.gz dbNSFP4.7a_variant.chr14.gz dbNSFP4.7a_variant.chr15.gz dbNSFP4.7a_variant.chr16.gz dbNSFP4.7a_variant.chr17.gz dbNSFP4.7a_variant.chr18.gz dbNSFP4.7a_variant.chr19.gz dbNSFP4.7a_variant.chr20.gz dbNSFP4.7a_variant.chr21.gz dbNSFP4.7a_variant.chr22.gz dbNSFP4.7a_variant.chrX.gz dbNSFP4.7a_variant.chrY.gz dbNSFP4.7a_variant.chrM.gz | grep -v "^#" ) | bgzip -c > dbNSFP4.7a_grch38.gz
```

### index the file with tabix
```
# version dbNSFP4.4a
tabix -s 1 -b 2 -e 2 dbNSFP4.4a_grch38.gz

# version dbNSFP4.7a (2024-05-16)
tabix -s 1 -b 2 -e 2 dbNSFP4.7a_grch38.gz
```

### remove the intermediate files
```
# version dbNSFP4.4a
rm dbNSFP4.4a_variant.chr*.gz
rm h dbNSFP4.4a.readme.txt search_dbNSFP47a.jar search_dbNSFP47a.readme.pdf search_dbNSFP47a.class tryhg19.in tryhg38.in LICENSE.txt tryhg18.in try.vcf

# version dbNSFP4.7a (2024-05-16)
rm dbNSFP4.7a_variant.chr*.gz
rm dbNSFP4.7a_header.txt
rm h dbNSFP4.7a.readme.txt search_dbNSFP47a.jar search_dbNSFP47a.readme.pdf search_dbNSFP47a.class tryhg19.in tryhg38.in LICENSE.txt tryhg18.in try.vcf
```

### for hg19
#### Replace coordinates by columns 8 and 9 (hg19 coordinates) and sort by those coordinates
#### TODO: explain the pipe command
#### or use something like here: https://github.com/GenomicsAotearoa/dbNSFP_build/blob/master/dbNSFP_pipeline_build.sh
```
# version dbNSFP4.4a
# transform hg38 to hg19
(zcat ../../hg38/dbnsfp/dbNSFP4.4a_grch38.gz | head -n 1 | sed -e 's/hg19/hg38/g'; zcat ../../hg38/dbnsfp/dbNSFP4.4a_grch38.gz | sed 1,1d | awk '{FS="\t"; OFS="\t"} $8 == "." { next; } $8 != "." { s=$1; $1=$8; $8=s; t=$2; $2=$9; $9=t; print; }' | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n) | bgzip -c > dbNSFP4.4a_grch37.gz

# Create tabix index
tabix -s 1 -b 2 -e 2 dbNSFP4.4a_grch37.gz

# version dbNSFP4.7a
# transform hg38 to hg19
(zcat ../../hg38/dbnsfp/dbNSFP4.7a_grch38.gz | head -n 1 | sed -e 's/hg19/hg38/g'; zcat ../../hg38/dbnsfp/dbNSFP4.7a_grch38.gz | sed 1,1d | awk '{FS="\t"; OFS="\t"} $8 == "." { next; } $8 != "." { s=$1; $1=$8; $8=s; t=$2; $2=$9; $9=t; print; }' | LC_ALL=C sort -t $'\t' -k1,1 -k2,2n) | bgzip -c > dbNSFP4.7a_grch37.gz

# Create tabix index
tabix -s 1 -b 2 -e 2 dbNSFP4.7a_grch37.gz
```
