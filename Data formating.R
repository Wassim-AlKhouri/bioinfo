library(maftools)        # to read MAF
library(GenomicRanges)   # for windows & counting
library(BSgenome.Hsapiens.UCSC.hg38)

######### Breast cancer SNV counts in 100 kb windows #########
maf <- read.maf('data/original/breast cancer (cohort)/cohortMAF.2025-04-28.maf')
unique_ids <- unique(maf@data$Tumor_Sample_Barcode)
head(unique_ids)

sample_id <- "TCGA-BH-A201-01A-11D-A14K-09" # Get our sample of interest

sample_maf <- subsetMaf(
  maf       = maf,
  tsb       = sample_id,       # tumor‐sample barcode
  mafObj    = TRUE             # return a maf object
)

snvs <- data.frame(
  chr = maf@data$Chromosome,
  pos = maf@data$Start_Position
)

# Create 100 kb tiling of the autosomes
bins <- tileGenome(
  seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg38)[paste0("chr",1:22)],  # to get only the autosomes
  tilewidth  = 100e3,
  cut.last.tile.in.chrom = TRUE # drop if smaller than 100 kb at the end of each chromo 
)

# Turn SNV list into a GRanges
snv_gr <- GRanges(
  seqnames = snvs$chr,
  ranges   = IRanges(start=snvs$pos, end = snvs$pos)
)

# Count SNVs per 100 kb bin
counts_vec <- countOverlaps(bins, snv_gr)

# Inspect the distribution
hist(log(counts_vec), breaks=30, main="SNVs per 100 kb window", xlab="Count")


df2 <- data.frame(
  seqnames   = as.character(seqnames(bins)),
  start      = start(bins),
  end        = end(bins),
  snv_count  = counts_vec
)

write.table(
  df2,
  file      =paste0("data/treated data/breast_cancer/snv_counts_100kb_windows_",sample_id,".tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

######### for arabidopsis #########

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

library(Biostrings)   # fast, memory‑efficient DNA manipulation
library(IRanges)      # simple range arithmetic
library(dplyr)        # nice tabular summaries

fasta_file <- 'data/original/Arabidopsis (TAIR)/TAIR10_chr_all.fas'
chr <- readDNAStringSet(fasta_file)
short_names        <- sub(" .*", "", names(chr))     # "Chr1", "Chr2", …
names(chr)    <- short_names

nuclear            <- chr[short_names %in% paste0("Chr", 1:5)]
nuclear

## count gc in 100 kb windows ##
gc_windows <- function(dna, chr_name, bin = 1e5) {
  starts  <- seq(1, length(dna), by = bin)
  ranges  <- IRanges(starts, width = pmin(bin, length(dna) - starts + 1))
  v       <- Views(dna, ranges)                       
  gc      <- rowSums(letterFrequency(v, c("G", "C")))
  
  tibble(chr       = chr_name,
         start     = start(ranges),
         end       = end(ranges),
         length_bp = width(ranges),
         GC_count  = gc,
         GC_prop    = gc/width(ranges) )
}

gc_table <- bind_rows(lapply(seq_along(nuclear), function(i)
  gc_windows(nuclear[[i]], names(nuclear)[i])))

print(gc_table, n = 8)

write.table(
  gc_table,
  file      =paste0('data/treated data/Arabidopsis/TAIR10_chr1-5_GC_100KB.tsv'),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE,
  col.names = TRUE
)