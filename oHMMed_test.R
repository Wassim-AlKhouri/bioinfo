################################################################################
# oHMMed on (i) breast-cancer SNV counts  |  (ii) Arabidopsis GC proportion   #
# --------------------------------------------------------------------------- #
#  paths assume you are inside the project root; tweak if needed              #
################################################################################

if (!requireNamespace("oHMMed", quietly = TRUE)) {
  install.packages("oHMMed")              # CRAN 2024-04-19
}
library(oHMMed)
library(ggmcmc)                           # plotting helpers

## ---------------------------------------------------------------------------
## 1.  LOAD THE WINDOW-LEVEL FEATURES YOU CREATED EARLIER
## ---------------------------------------------------------------------------
snv_df <- read.table(
  "data/treated data/breast_cancer/snv_counts_100kb_windows_TCGA-BH-A201-01A-11D-A14K-09.tsv",
  header = TRUE, sep = "\t"
)
counts_raw <- as.integer(snv_df$snv_count)        # make sure they are *integers*
counts_raw <- pmax(counts_raw, 1L)           # avoid zero counts (log(0) = -Inf)

gc_df <- read.table(
  "data/treated data/Arabidopsis/TAIR10_chr1-5_GC_100KB.tsv",
  header = TRUE, sep = "\t"
)

# ---- NEW LINE ----  unify column header with SNV table
if ("chr" %in% names(gc_df)) names(gc_df)[names(gc_df) == "chr"] <- "seqnames"

gc_prop <- gc_df$GC_prop


## ---------------------------------------------------------------------------
## 2.  CHOOSE NUMBER OF HIDDEN STATES  +  BUILD REASONABLE PRIORS
## ---------------------------------------------------------------------------
K_counts <- 4                     # try 4 for mutation-rate landscape
K_gc     <- 3                     # try 3 for GC (low-mid-high)

### 2a. Priors for the **gamma-Poisson** model  (counts)
prior_T_counts   <- generate_random_T(K_counts)    # near-diagonal random matrix
prior_alpha      <- 3                              # single shape prior (fairly vague)

# pick K quantiles of the observed counts as prior means,
# then convert means to (descending) rate parameters  beta = alpha / mean
mean_q           <- as.numeric(quantile(counts_raw,
                               probs = seq(0.2, 0.8, length.out = K_counts)))
prior_betas      <- prior_alpha / sort(mean_q, decreasing = TRUE)

### 2b. Priors for the **normal** model  (GC)
prior_T_gc       <- generate_random_T(K_gc)
prior_means      <- as.numeric(quantile(gc_prop,
                               probs = seq(0.2, 0.8, length.out = K_gc)))
prior_sd         <- sd(gc_prop)                    # single shared σ prior


## ---------------------------------------------------------------------------
## 3.  RUN oHMMed
## ---------------------------------------------------------------------------
set.seed(42)
fit_snv <- hmm_mcmc_gamma_poisson(
  data         = counts_raw,
  prior_T      = prior_T_counts,
  prior_betas  = prior_betas,
  prior_alpha  = prior_alpha,
  iter         = 80,
  warmup       = 20,
  thin         = 10,
  print_params = FALSE,     # turn on if you like to watch the chain
  verbose      = TRUE
)

set.seed(43)
fit_gc <- hmm_mcmc_normal(
  data         = gc_prop,
  prior_T      = prior_T_gc,
  prior_means  = prior_means,
  prior_sd     = prior_sd,
  iter         = 80,
  warmup       = 20,
  thin         = 10,
  print_params = FALSE,
  verbose      = TRUE
)


## ---------------------------------------------------------------------------
## 4.  QUICK DIAGNOSTICS  (optional but highly recommended)
## ---------------------------------------------------------------------------
plot(fit_snv)                                # trace- & density plots, LL trace, etc.
plot(fit_gc)

# If you prefer ggmcmc style:
ggmcmc::ggs_traceplot(convert_to_ggmcmc(fit_snv, pattern = "beta"))
ggmcmc::ggs_traceplot(convert_to_ggmcmc(fit_gc,  pattern = "mean"))


##############################################################################
# 5.  MOST-LIKELY STATE PER WINDOW  +  WRITE BED-LIKE TRACKS
##############################################################################

## ---------- SNV counts (gamma-Poisson) ------------------------------------
state_snv <- fit_snv$estimates$posterior_states               # factor (L)
post_snv  <- fit_snv$estimates$posterior_states_prob          # L × K matrix

snv_df$state <- state_snv

write.table(
  snv_df[, c("seqnames", "start", "end", "snv_count", "state")],
  file      = "data/treated data/breast_cancer/TCGA-BH-A201-01A-11D-A14K-09_oHMMed_K4.bed",
  sep       = "\t", quote = FALSE, row.names = FALSE
)

## ---------- Arabidopsis GC (normal) ---------------------------------------
state_gc <- fit_gc$estimates$posterior_states
post_gc  <- fit_gc$estimates$posterior_states_prob

gc_df$state <- state_gc

write.table(
  gc_df[, c("seqnames", "start", "end", "GC_prop", "state")],
  file      = "data/treated data/Arabidopsis/TAIR10_chr1-5_GC_100KB_oHMMed_K3.bed",
  sep       = "\t", quote = FALSE, row.names = FALSE
)



## ---------------------------------------------------------------------------
## 6.  OPTIONAL: FIT ALTERNATIVE K AND COMPARE
## ---------------------------------------------------------------------------
#altK <- lapply(c(3,5), function(k) {
#  hmm_mcmc_gamma_poisson(counts_raw,
#                         prior_T = generate_random_T(k),
#                         prior_betas = rep(mean(prior_betas), k),
#                         prior_alpha = prior_alpha,
#                         iter = 80, warmup = 20, thin = 10,
#                         print_params = FALSE, verbose = FALSE)
#})
#K3 <- altK[[1]] ; K5 <- altK[[2]]
#
#comp_snv <- data.frame(
#  K        = c(3,4,5),
#  logLik   = c(mean(K3$estimates$log_likelihood),
#               mean(fit_snv$estimates$log_likelihood),
#               mean(K5$estimates$log_likelihood)),
#  entropy  = c(K3$estimates$segmentation_entropy,
#               fit_snv$estimates$segmentation_entropy,
#               K5$estimates$segmentation_entropy)
#)
#print(comp_snv, row.names = FALSE)