



################# Article verification ############

df <- read.table("~/Desktop/bioinf/M2/methods bionfo /GenomeAnnotations/MouseAnnotation.txt",
                 header = TRUE, 
                 sep="")
gc <- df$GC_PROPORTION
gc <- na.omit(gc)
var(gc)
d_gc  <- diff(gc)           # ΔGC = GC[i+1] − GC[i]

# distribution 

hist(gc,
     breaks = 100,               # number of bins; adjust to taste
     freq   = FALSE,           # density rather than counts
     col    = "lightgray",
     main   = "Distribution of SNV-count (Arabidopsis) ",
     xlab   = "log SNV count")







# expecting same result has in the article mouse var = 0.00055 
var_d_gc  <- var(d_gc)
var_d_gc
set.seed(123)
d_rand   <- diff(sample(gc))
var_rand <- var(d_rand)
var_rand
ft <- var.test(d_gc, d_rand)
print(ft)

############ Breast cancer #############

data <- read.table('/Volumes/T7 Shield/oHMMed /DATA/BREAST/100KB_data /snv_counts_100kb_windows_TCGA-BH-A201-01A-11D-A14K-09.tsv',
                   sep = "\t",header = TRUE )
counts <- data$snv_count
counts_log <- log1p(counts)



# distribution

hist(counts_log,
     breaks = 20,               # number of bins; adjust to taste
     freq   = FALSE,           # density rather than counts
     col    = "lightgray",
     main   = "Distribution of SNV-count ",
     xlab   = "log SNV count")


d_obs   <- diff(counts)
var_obs <- var(d_obs)







for(i in 1:5){
  set.seed(123+i)
  d_rand   <- diff(sample(counts))
  var_rand <- var(d_rand)
  ft <- var.test(d_obs, d_rand)
  print(ft)
}

# we want var_obs < var_rand  and p-val << 0.05 -> 
# alternative hypothesis: true ratio of variances is not equal to 1




####################### arabidopsis ##################



data <- read.table('/Users/godinmax/Desktop/bioinf/M2/methods bionfo /MAX_WAS/bioinfo/TAIR10_chr1-5_GC_100KB.tsv',
                   sep = "\t",header = TRUE )
counts <- data$GC_prop

hist(counts,
     breaks = 100,               # number of bins; adjust to taste
     freq   = FALSE,           # density rather than counts
     col    = "lightgray",
     main   = "Distribution of SNV-count (Arabidopsis) ",
     xlab   = "SNV count")


d_obs   <- diff(counts)
var_obs <- var(d_obs)
var_obs

for(i in 1:5){
  set.seed(123+i)
  d_rand   <- diff(sample(counts))
  var_rand <- var(d_rand)
  ft <- var.test(d_obs, d_rand)
  print(ft)
}
