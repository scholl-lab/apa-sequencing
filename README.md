# Aldosterone-Producing Adenoma sequencing project (APA-sequencing)

Code for the aldosterone-producing adenoma (APA) sequencing project.
It also contains general descriptions about the project, e.g. aim, sample collection and library preparation.

## Documentation

Further project documentation on [GitHub Pages](https://scholl-lab.github.io/apa-sequencing/).

## List of currently known APA-associated genes

- KCNJ5
- CACNA1D
- ATP1A1
- ATP2B3
- CLCN2
- CACNA1H
- SLC30A1
- CTNNB1
- GNAQ
- GNA11
- CADM1
- GNAS
- MCOLN3

## TODO other

- [ ] add information about Sanger sequencing of known APA-associated genes
- [ ] list off all samples included in the project including the ones screened out through Sanger sequencing of KCNJ5

## TODO analyses

- SNV/indel variants in known APA-associated genes to find causative mutations in these genes
- SNV/indel variants in other genes to find novel APA-associated genes
  1. filter COSMIC for samples with adenoma producing tumors and look for recurrent mutations or frequently mutated genes
  2. look for rare disease associations with the HPO terms Neoplasm of the adrenal gland HP:0100631 (or more specific: Neoplasm of the adrenal cortex HP:0100641) or Hyperaldosteronism HP:0000859 (or more specific: Primary hyperaldosteronism HP:0011736)
  3. interaction analysis of APA-associated genes with other genes using STRING
  4. exome wide analysis of somatic SNV/indel variants and comparison with other cohorts of APA tumors
  5. mutational rate for these tumors
- CNV analysis
  1. recurrently affected regions in APA tumors
  2. comparison with other cohorts of APA tumors (review published data)
- accurate calling (freebayes, other callers) for known APA-associated driver mutations (Example Sanger: KCNJ5)
- genomes
  1. SNV/indel and CNV like for exomes
  2. comparison with exome data
- non-functional (NF) vs. functional APA tumors (F), these samples have genome and exome data
  1. pairwise calling of SNV/indel variants in functional vs. non-functional APA tumors (F vs. NF, NF vs. F, F vs. N, NF vs. N) to find differences in the mutational landscape especially secondary hits that cause proliferation or hormone production
  2. evolutionary analysis of functional vs. non-functional APA tumors to find the order of mutations


## filter with variantcentrifuge
variantcentrifuge --reference GRCh38.p14 -G testing/apa_and_cosmic_gene_list.txt -v testing/APA79/APA79-T_To.filtered.annotated.vcf.gz --log-level DEBUG --xlsx --html-report -f "((((dbNSFP_gnomAD_exomes_AC[0] <= 2 ) | (na dbNSFP_gnomAD_exomes_AC[0])) & (( dbNSFP_gnomAD_genomes_AC[0] <= 2 ) | (na dbNSFP_gnomAD_genomes_AC[0]))) | ((exists ID) & ( ID =~ 'COS' ))) & (( na FILTER ) | (FILTER = 'PASS')) & ( GEN[*].DP >= 50 ) & ( GEN[*].AF >= 0.05 ) & ( GEN[*].AF < 0.35 ) & ( VARTYPE = 'SNP' )" --preset high_or_moderate --threads 16 --append-extra-sample-fields GEN[*].DP GEN[*].AD GEN[*].AF --igv --igv-reference hg38_1kg --bam-mapping-file testing/apa_bam_mapping.txt

variantcentrifuge --reference GRCh38.p14 -g all -v testing/APA79/APA79_TvsN.filtered.annotated.vcf.gz --log-level DEBUG --xlsx --html-report -f "((((dbNSFP_gnomAD_exomes_AC[0] <= 2 ) | (na dbNSFP_gnomAD_exomes_AC[0])) & (( dbNSFP_gnomAD_genomes_AC[0] <= 2 ) | (na dbNSFP_gnomAD_genomes_AC[0]))) | ((exists ID) & ( ID =~ 'COS' ))) & (( na FILTER ) | (FILTER = 'PASS')) & ( GEN[*].DP >= 50 ) & ( GEN[0].AF < 0.03 ) & ( GEN[1].AF >= 0.05 ) & ( GEN[1].AF < 0.45 ) & ( VARTYPE = 'SNP' )" --preset high_or_moderate --threads 16 --append-extra-sample-fields GEN[*].DP GEN[*].AD GEN[*].AF --igv --igv-reference hg38_1kg --bam-mapping-file testing/apa_bam_mapping.txt


## qualimap
qualimap bamqc --java-mem-size=8G -c -gff /data/cephfs-1/work/groups/scholl/shared/target_files/agilent/S33266340/hg38/S33266340_Regions_noheader.bed -bam results/exomes/bqsr/APA79-N.merged.dedup.bqsr.bam


## subset BAMs
ls results/exomes/bqsr/*.bam | parallel -j +0 "samtools view -F 1024 {} -L /data/cephfs-1/work/groups/scholl/shared/target_files/apa_genes/hg38/apa_genes.genes2bed.GRCh38.S33266436_Regions.padding1000bp.bed -bo results/exomes/samtools_view_apa_genes_padding/{/.}.apa-genes.bam"

ls results/exomes/samtools_view_apa_genes_padding/*.bam | parallel -j +0 "samtools index {}"

## indexcov
goleft indexcov --directory results/exomes/qc/indexcov results/exomes/bqsr/*.bam