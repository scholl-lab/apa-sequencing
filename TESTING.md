
zcat download/lane36/P1775_DNA_03_iAPA1-N_S3_L002_R1_001.fastq.gz | head -n 50000 | gzip > download/P1775_DNA_03_iAPA1-N_S3_L002_R1_001.head50k.fastq.gz 
zcat download/lane36/P1775_DNA_03_iAPA1-N_S3_L002_R2_001.fastq.gz | head -n 50000 | gzip > download/P1775_DNA_03_iAPA1-N_S3_L002_R2_001.head50k.fastq.gz


bwa mem -t 10 -R "@RG\tID:lane36-S3-L002-iAPA1-N\tSM:iAPA1-N\tLB:iAPA1-N\tPL:ILLUMINA\tPU:lane36-S3-L002" analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz download/P1775_DNA_03_iAPA1-N_S3_L002_R1_001.head50k.fastq.gz download/P1775_DNA_03_iAPA1-N_S3_L002_R2_001.head50k.fastq.gz | samtools sort -O BAM -o results/P1775_DNA_03_iAPA1-N_S3_L002_lane36.head50k.bam -



zcat download/lane36/P1775_DNA_09_iAPA3-N_S9_L002_R1_001.fastq.gz | head -n 5000000 | gzip > download/P1775_DNA_09_iAPA3-N_S9_L002_R1_001.head5000k.fastq.gz
zcat download/lane36/P1775_DNA_09_iAPA3-N_S9_L002_R2_001.fastq.gz | head -n 5000000 | gzip > download/P1775_DNA_09_iAPA3-N_S9_L002_R2_001.head5000k.fastq.gz

bwa mem -t 12 -R "@RG\tID:lane36-S9-L002-iAPA3-N\tSM:iAPA3-N\tLB:iAPA3-N\tPL:ILLUMINA\tPU:lane36-S9-L002" analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz download/P1775_DNA_09_iAPA3-N_S9_L002_R1_001.head500k.fastq.gz download/P1775_DNA_09_iAPA3-N_S9_L002_R2_001.head500k.fastq.gz | samtools sort -@4 -m4g -O BAM -o results/aligned/P1775_DNA_09_iAPA3-N_S9_L002_lane36.head500k.bam - 


# MarkDuplicatesSpark with gatk
# based on:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
# https://gatk.broadinstitute.org/hc/en-us/articles/13832682540699-MarkDuplicatesSpark
gatk MarkDuplicatesSpark \
    -I results/merged/P1775_DNA_03_iAPA1-N.merged.bam \
    -O P1775_DNA_03_iAPA1-N.merged.dedup.bam \
    -M P1775_DNA_03_iAPA1-N.merged.dedup_metrics.txt \
    --conf 'spark.executor.cores=8'

# BaseRecalibrator with gatk
# based on:
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery
# https://gatk.broadinstitute.org/hc/en-us/articles/13832694104987-BaseRecalibratorSpark-BETA-
# https://gatk.broadinstitute.org/hc/en-us/articles/13832711710107-ApplyBQSRSpark-BETA-
gatk BaseRecalibratorSpark \
    -I results/dedup/P1775_DNA_01_iAPA1-F.merged.dedup.bam \
    -R analysis/ref/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    --known-sites analysis/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites analysis/GATK_resource_bundle/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites analysis/GATK_resource_bundle/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    -O results/bqsr/P1775_DNA_01_iAPA1-F.merged.dedup.recal_data.table \
    --conf 'spark.executor.cores=8'

gatk ApplyBQSRSpark \
    -I results/dedup/P1775_DNA_01_iAPA1-F.merged.dedup.bam \
    -bqsr results/bqsr/P1775_DNA_01_iAPA1-F.merged.dedup.recal_data.table \
    -O results/bqsr/P1775_DNA_01_iAPA1-F.merged.dedup.bqsr.bam \
    --conf 'spark.executor.cores=8'

# check out snakemake wrappers
#
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/gatk/markduplicatesspark.html

## GATK with conda in snakemake
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html


# mutect2 
# have to parallelize into intervalls somehow, will be slow else
# --> PON: use GATK bundle
# --> intervals: use GATK bundle if matching reference, otherwise make own intervals
# quick solution: make an unparallelized Mutect script and start it already (in case the parallelized one takes long to write)
https://github.com/berntpopp/conNDD/blob/main/tertiary_analysis/Commands_Alignment-And-Calling.sh
https://gatk.broadinstitute.org/hc/en-us/community/posts/360060460131-GATK4-Parallelizing-genotypegvcfs
https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists
https://gatk.broadinstitute.org/hc/en-us/community/posts/360062233271-Parallelizing-Mutect2
https://gatk.broadinstitute.org/hc/en-us/articles/360035531132
https://gatk.broadinstitute.org/hc/en-us/articles/13832710384155-Mutect2


# CNV GATK
https://gatk.broadinstitute.org/hc/en-us/articles/360035531092#2

# cnv CNVkit
https://cnvkit.readthedocs.io/en/stable/nonhybrid.html
--> output as vcf possible?

# annotation
https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial
