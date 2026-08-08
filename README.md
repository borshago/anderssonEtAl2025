This repository contains the RNA-seq data processing and analysis files for the publication Andersson et al. (2025). Old mitochondria regulate niche renewal via α-ketoglutarate metabolism in stem cells. Nat Metab 7, 1344–1357. https://doi.org/10.1038/s42255-025-01325-7.

Contents:

	- environment.yml: specification for creating the conda environment
	- environment_export.yml: the full exported conda environment
	- Snakefile: snakemake workflow for data processing and analysis
	- rename.sh: bash script for renaming the original FASTQ files
	- de_gsea_analysis.R: script for differential expression and gene set enrichment analyses
	- de_gsea_panethContam.R: script for differential expression and gene set enrichment analyses comparing ISCmito-Ys and ISCmito-Ys with deliberate 5% Paneth cell contamination


The scripts expect the following directory structure:

```
  CODE/
  DATA/
  ├─ ALIGNED/
  ├─ COUNTS/
  ├─ GENOME/
  │  ├─ STAR/
  ├─ PROCESSED/
  RESULTS/
  ├─ OLD_VS_YOUNG/
  │  ├─ DE/
  │  ├─ GSEA/
  ├─ PANETH_CONTAM/
  │  ├─ DE/
  │  ├─ GSEA/
  ├─ QCOUT/
  │  ├─ RSEQC/
  environment.yml
  Snakefile
```
