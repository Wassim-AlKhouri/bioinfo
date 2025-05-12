



################# Article verification ############

df <- read.table("data/original/GenomeAnnotations (paper)/GenomeAnnotations/MouseAnnotation.txt",
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
     main   = "Distribution of SNV-count (mouse) ",
     xlab   = "GC count")

set.seed(2)
gc.sub <- sample(gc, 5000)
shapiro.test(gc.sub)




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

data <- read.table('data/treated data/breast_cancer/snv_counts_100kb_windows_TCGA-BH-A201-01A-11D-A14K-09.tsv',
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

hist(d_obs,
     breaks = 100,               # number of bins; adjust to taste
     freq   = FALSE,           # density rather than counts
     col    = "lightgray",
     main   = "Pairwise differences in average SNV-count ",
     xlab   = "differences in SNV-count frequency ")

var_obs <- var(d_obs)
var_obs

set.seed(123)
d_rand   <- diff(sample(counts))
var_rand <- var(d_rand)
ft <- var.test(d_obs, d_rand)
print(var_rand)
print(ft)
hist(d_rand,
     breaks = 50,
     freq   = FALSE,
     col    = rgb(1,1,1,0.4),   # semi-transparent white
     border = "white",          # white outline
     add    = TRUE)



for(i in 1:5){
  set.seed(123+i)
  d_rand   <- diff(sample(counts))
  var_rand <- var(d_rand)
  ft <- var.test(d_obs, d_rand)
  print(var_rand)
  print(ft)
}

# we want var_obs < var_rand  and p-val << 0.05 -> 
# alternative hypothesis: true ratio of variances is not equal to 1




####################### arabidopsis ##################



data <- read.table('data/treated data/Arabidopsis/TAIR10_chr1-5_GC_100KB.tsv',
                   sep = "\t",header = TRUE )
counts <- data$GC_prop

hist(counts,
     breaks = 50,               # number of bins; adjust to taste
     freq   = FALSE,           # density rather than counts
     col    = "lightgray",
     main   = "Distribution of SNV-count ",
     xlab   = "GC count")


d_obs   <- diff(counts)

var_obs <- var(d_obs)
var_obs

hist(d_obs,
     breaks = 50,               # number of bins; adjust to taste
     freq   = FALSE,           # density rather than counts
     col    = "lightgray",
     main   = "Pairwise differences in average GC proportion  ",
     xlab   = "differences in GC proportion ")


set.seed(123)
d_rand   <- diff(sample(counts))
var_rand <- var(d_rand)
print(var_rand)
ft <- var.test(d_obs, d_rand)
print(ft)

hist(d_rand,
     breaks = 50,
     freq   = FALSE,
     col    = rgb(1,1,1,0.4),   # semi-transparent white
     border = "white",          # white outline
     add    = TRUE)



for(i in 1:5){
  set.seed(123+i)
  d_rand   <- diff(sample(counts))
  var_rand <- var(d_rand)
  print(var_rand)
  ft <- var.test(d_obs, d_rand)
  print(ft)
}
