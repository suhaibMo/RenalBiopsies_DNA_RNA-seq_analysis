# RenalBiopsies_DNA_RNA-seq_analysis
Integrated multi-omics analysis for kidney transplant patients to decipher molecular gene signatures for delayed graft function. 

**Repository sructure**

## RNA-seq analysis
- This repository contains script thats reads all raw illumina fastq files from a folder and performs quality control, genome alignment and read summarisation to quantify read counts
- Genome alignment and read summarisation (gene quantification)
- Trim galore for trimming low quality reads fastq files
- Tophat/Bowtie for aligning reads to the updated reference genome (GRCh38) from ensemble 
- Scripts for analysis of RNA-seq data for differential gene expression (DGE) using DESeq2

## BS-seq analysis
Scripts for analysis of DNA-seq to extract methylation changes (epigenetic) footprints on CpG islands 

## Epigenetics
Scripts for mapping differential gene expression signatures and associated methylation status

## Glossary notation on samples
- b - pre perfused kidney patients
- b1 - post perfused kidney patients
- DGF - Delayed graft function conditioned patients
- IGF - No Delayed graft function patients 

The analysis of these results are published in Aging cell [A molecular signature for delayed graft function](https://onlinelibrary.wiley.com/doi/full/10.1111/acel.12825 "paper")

```
McGuinness, D., Mohammed, S., Monaghan, L., Wilson, P. A., Kingsmore, D. B., Shapter, O., Stevenson, K. S., Coley, S. M., Devey, L., Kirkpatrick, R. B., & Shiels, P. G. (2018). A molecular signature for delayed graft function. Aging Cell, 17(5), e12825. https://doi.org/10.1111/acel.12825
```