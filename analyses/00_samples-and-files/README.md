# General project overview including samples and files

## Sample Metadata and Known Issues

### Known Sample Swaps
The project includes several known sample swaps that are corrected in `calling-metadata.R`:

1. **APA12/APA13 Tumor Swap**: The tumor samples for patients APA12 and APA13 were swapped. 
   - The script swaps sample names and BAM file basenames while maintaining individual identifiers.

2. **APA58 T/N Swap**: The tumor and normal samples for patient APA58 were swapped.
   - Sample designations are swapped and analysis type is standardized to TvsN.

3. **iAPA5/iAPA6 Complex Reciprocal Swap**: A complex swap between iAPA5 and iAPA6 samples:
   - iAPA5 group contains: iAPA6-F, iAPA6-N, and iAPA5-NF samples
   - iAPA6 group contains: iAPA5-F, iAPA5-N, and iAPA6-NF samples
   - For each individual, the F and N samples have swapped prefixes while NF samples remain with their original prefix.

### Sample Exclusions
Several samples are excluded from analysis due to quality issues, as detailed in the `calling-metadata.R` script.

## TODO tickets
- [ ] metadata for genomes project part
- [ ] unify metadata for exomes and genomes project parts
- [ ] add md5sums for all files from sequencing provider and after download