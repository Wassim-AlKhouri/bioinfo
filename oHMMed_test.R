###############################################################################
# oHMMed example workflow
#  ├─ (i) SNV counts in 100-kb breast-cancer windows  → gamma–Poisson HMM
#  └─ (ii) Arabidopsis GC proportion                → Gaussian HMM
# -----------------------------------------------------------------------------
#  The script is robust to K ≥ 2 … 9; empty-state crashes are avoided by
#  (1) using a strictly stochastic, near-diagonal transition matrix and
#  (2) seeding new states with sensible emission priors.
###############################################################################

## ────────────────────────────────────────────────────────────────────────────
## 0.  PACKAGES
## ────────────────────────────────────────────────────────────────────────────
req_pkgs <- c("oHMMed", "dplyr", "tidyr", "ggplot2", "ggmcmc")
for (pkg in req_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

## ────────────────────────────────────────────────────────────────────────────
## 1.  DATA
## ────────────────────────────────────────────────────────────────────────────
# ---- (i) SNV counts ────────────────────────────────────────────────────────
snv_df <- read.table(
  "data/treated data/breast_cancer/snv_counts_100kb_windows_TCGA-BH-A201-01A-11D-A14K-09.tsv",
  header = TRUE, sep = "\t"
)
counts_raw <- as.integer(snv_df$snv_count)
counts_raw <- pmax(counts_raw, 1L)               # avoid log(0) = –Inf

# ---- (ii) Arabidopsis GC proportion ────────────────────────────────────────
gc_df <- read.table(
  "data/treated data/Arabidopsis/TAIR10_chr1-5_GC_100KB.tsv",
  header = TRUE, sep = "\t"
)
if ("chr" %in% names(gc_df)) names(gc_df)[names(gc_df) == "chr"] <- "seqnames"
gc_prop <- gc_df$GC_prop

## ────────────────────────────────────────────────────────────────────────────
## 2.  HELPER: SAFE RANDOM TRANSITION MATRIX
## ────────────────────────────────────────────────────────────────────────────
generate_random_T <- function(k, self = 0.85) {
  stopifnot(k >= 2, self > 0, self < 1)

  T <- matrix(0, k, k)
  remain <- 1 - self

  for (i in 1:k) {
    T[i, i] <- self
    if (i > 1)  T[i, i - 1] <- remain * ifelse(i == k, 1, 0.5)
    if (i < k)  T[i, i + 1] <- remain * ifelse(i == 1, 1, 0.5)
  }

  # (optional) sanity check – every row must be 1
  if (any(abs(rowSums(T) - 1) > 1e-12))
    stop("Internal error: row sums are not 1")

  T
}

## ────────────────────────────────────────────────────────────────────────────
## 3.  CHOOSE K AND BUILD PRIORS
## ────────────────────────────────────────────────────────────────────────────
K_counts <- 4               # try 4 first; increase only after diagnostics
K_gc     <- 4

# --- gamma–Poisson priors (SNV counts) --------------------------------------
prior_alpha      <- 3
mean_q           <- quantile(counts_raw,
                             probs = seq(0.1, 0.9, length.out = K_counts),
                             names = FALSE)
prior_betas      <- prior_alpha / sort(pmax(mean_q, 10), decreasing = TRUE)
prior_T_counts   <- generate_random_T(K_counts, self = 0.85)

# --- Gaussian priors (GC proportion) ----------------------------------------
prior_sd         <- sd(gc_prop)
prior_means      <- sort(quantile(gc_prop,
                                  probs = seq(0.1, 0.9, length.out = K_gc),
                                  names = FALSE))
prior_T_gc       <- generate_random_T(K_gc, self = 0.85)

## ────────────────────────────────────────────────────────────────────────────
## 4.  RUN oHMMed
## ────────────────────────────────────────────────────────────────────────────
options(error = recover)    # drop into debugger on any unexpected error

set.seed(42)
fit_snv <- hmm_mcmc_gamma_poisson(
  data         = counts_raw,
  prior_T      = prior_T_counts,
  prior_betas  = prior_betas,
  prior_alpha  = prior_alpha,
  iter         = 100,
  warmup       =   20,
  thin         =    10,
  print_params = FALSE,
  verbose      = TRUE
)

set.seed(43)
fit_gc <- hmm_mcmc_normal(
  data         = gc_prop,
  prior_T      = prior_T_gc,
  prior_means  = prior_means,
  prior_sd     = prior_sd,
  iter         = 100,
  warmup       =   20,
  thin         =    10,
  print_params = FALSE,
  verbose      = TRUE
)

##############################################################################
# 5.  DIAGNOSTICS  +  SAVE EVERY FIGURE
##############################################################################
output_dir <- "data/treated data/plots"          # <- where images will go
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ---------- 5a. built-in diagnostic plots ----------------------------------
# oHMMed::plot() returns a grid object; wrap it in a PNG device
png(file.path(output_dir, "diagnostics_snv_K4.png"),
    width = 1600, height = 1200, res = 150)
plot(fit_snv, main = "oHMMed diagnostics — SNV (K = 4)")
dev.off()

png(file.path(output_dir, "diagnostics_gc_K4.png"),
    width = 1600, height = 1200, res = 150)
plot(fit_gc,  main = "oHMMed diagnostics — GC (K = 4)")
dev.off()

## ---------- 5b. ggmcmc trace plots -----------------------------------------
library(ggmcmc)

p_beta  <- ggs_traceplot(convert_to_ggmcmc(fit_snv, pattern = "beta"))
p_mean  <- ggs_traceplot(convert_to_ggmcmc(fit_gc,  pattern = "mean"))

ggplot2::ggsave(filename = file.path(output_dir, "trace_beta_snv_K4.png"),
                plot = p_beta, width = 9, height = 6, dpi = 150)
ggplot2::ggsave(filename = file.path(output_dir, "trace_mean_gc_K4.png"),
                plot = p_mean, width = 9, height = 6, dpi = 150)

## ────────────────────────────────────────────────────────────────────────────
## 6.  MOST-LIKELY STATE PER WINDOW  →  BED-LIKE TRACKS
## ────────────────────────────────────────────────────────────────────────────
snv_df$state <- fit_snv$estimates$posterior_states
write.table(
  snv_df[, c("seqnames", "start", "end", "snv_count", "state")],
  file      = "data/treated data/breast_cancer/TCGA-BH-A201-01A-11D-A14K-09_oHMMed_K4.bed",
  sep       = "\t", quote = FALSE, row.names = FALSE
)

gc_df$state <- fit_gc$estimates$posterior_states
write.table(
  gc_df[, c("seqnames", "start", "end", "GC_prop", "state")],
  file      = "data/treated data/Arabidopsis/TAIR10_chr1-5_GC_100KB_oHMMed_K4.bed",
  sep       = "\t", quote = FALSE, row.names = FALSE
)

##############################################################################
# 7.  SWEEP OVER K   +   SAVE EVERY EXTRA PLOT
##############################################################################
states <- 2:9
fits_gc <- vector("list", length(states))

for (i in seq_along(states)) {
  k <- states[i]
  message("◆ Fitting K = ", k)
  set.seed(100 + k)
  fits_gc[[i]] <- hmm_mcmc_normal(
    data         = gc_prop,
    prior_T      = generate_random_T(k, self = 0.85),
    prior_means  = sort(quantile(gc_prop,
                                 probs = seq(0.1, 0.9, length.out = k),
                                 names = FALSE)),
    prior_sd     = prior_sd,
    iter         = 800, warmup = 160, thin = 10,
    print_params = FALSE, verbose = FALSE
  )

  ## save the diagnostic plot for this K
  png(file.path(output_dir,
                sprintf("diagnostics_gc_K%02d.png", k)),
      width = 1600, height = 1200, res = 150)
  plot(fits_gc[[i]],
       main = sprintf("oHMMed diagnostics — GC, K = %d", k))
  dev.off()
}

## --- log-likelihood summary figure -----------------------------------------
loglik_gc_summary <- data.frame(
  K        = states,
  MeanLL   = vapply(fits_gc, function(f) mean(f$estimates$log_likelihood),
                    numeric(1)),
  MedianLL = vapply(fits_gc, function(f) median(f$estimates$log_likelihood),
                    numeric(1))
)

library(ggplot2)
p_ll <- ggplot(loglik_gc_summary, aes(x = K)) +
          geom_line(aes(y = MeanLL),    linewidth = 1) +        # <- linewidth!
          geom_point(aes(y = MeanLL),   size = 2) +
          geom_line(aes(y = MedianLL),  linetype = "dashed", linewidth = 1) +
          geom_point(aes(y = MedianLL), size = 2) +
          scale_x_continuous(breaks = states) +
          labs(
            title = "Arabidopsis GC: mean / median log-likelihood vs. K",
            x     = "Number of hidden states (K)",
            y     = "Log-likelihood"
          ) +
          theme_bw()

ggsave(filename = file.path(output_dir, "gc_loglik_vs_K.png"),
       plot = p_ll, width = 9, height = 6, dpi = 150)