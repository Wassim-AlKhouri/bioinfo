# oHMMed Genomic Segmentation Pipeline

## Purpose  
Segment genomic windows into ordered hidden states with **oHMMed**.  
Outputs:

* 100 kb tracks of SNV burden and GC content  
* Diagnostic convergence plots  
* BED files ready for genome browsers  

---

## Repository layout

| Path | Content |
|------|---------|
| `oHMMed.R` | End-to-end model fitting + plots |
| `oHMMed assumptions verif.R` | Distribution checks vs. paper |
| `Data formating.R` | Raw FASTA/MAF → 100 kb tables |
| `data/original/` | Raw inputs (FASTA, MAF, annotations) |
| `data/treated data/` | Derived tables, plots, BED tracks |
| `README.md` | **← this file** |

---

## Requirements
### R packages

| Source | Packages |
|--------|----------|
| CRAN | `oHMMed`, `dplyr`, `tidyr`, `ggplot2`, `ggmcmc`, `maftools` |
| Bioconductor | `Biostrings`, `IRanges`, `GenomicRanges`,<br>`BSgenome.Hsapiens.UCSC.hg38`, `BSgenome.Athaliana.TAIR.TAIR10` |

---

## Quick install

```r
# Bioconductor (one-off)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN
install.packages(c(
  "oHMMed", "dplyr", "tidyr", "ggplot2", "ggmcmc", "maftools"
))

# Bioconductor
BiocManager::install(c(
  "Biostrings", "IRanges", "GenomicRanges",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Athaliana.TAIR.TAIR10"
))
```

## Input data

| Script | Needed files |
|--------|--------------|
| `Data formating.R` | *Breast cancer* MAF → `data/original/breast cancer (cohort)/cohortMAF*.maf`<br>*Arabidopsis* FASTA → `data/original/Arabidopsis (TAIR)/TAIR10_chr_all.fas` |
| `oHMMed assumptions verif.R` | Mouse annotation → `data/original/GenomeAnnotations (paper)/GenomeAnnotations/MouseAnnotation.txt` |
| `oHMMed.R` | Uses the tables created by `Data formating.R` |

Decompress the data in `data/original/breast cancer (cohort)` and `data/original/Arabidopsis (TAIR)` directories before running the scripts.