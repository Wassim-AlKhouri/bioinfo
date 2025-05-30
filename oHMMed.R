## ────────────────────────────────────────────────────────────────────────────
## How to run.
## 1.  Install the required packages (if not already installed).
## 2.  Extract the compressed data files in the 'data/original' directory (in their respective subdirectories).
## 3.  Run this script in R or RStudio.
## ────────────────────────────────────────────────────────────────────────────

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
# ───── (i) SNV counts ────────────────────────────────────────────────────────
snv_df <- read.table(
  "data/treated data/breast_cancer/snv_counts_100kb_windows_TCGA-BH-A201-01A-11D-A14K-09.tsv",
  header = TRUE, sep = "\t"
)
counts_raw <- as.integer(snv_df$snv_count)
counts_raw <- pmax(counts_raw, 1L)               # avoid log(0) = –Inf

# ───── (ii) Arabidopsis GC proportion ────────────────────────────────────────
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

  # sanity check – every row must be 1
  if (any(abs(rowSums(T) - 1) > 1e-12))
    stop("Internal error: row sums are not 1")

  T
}

## ────────────────────────────────────────────────────────────────────────────
## 3.  CHOOSE K AND BUILD PRIORS
## ────────────────────────────────────────────────────────────────────────────
K_counts <- 4               # k for test run (not final K)
K_gc     <- 4

# ───── gamma–Poisson priors (SNV counts) ─────────────────────────────────────
prior_alpha      <- 3
mean_q           <- quantile(counts_raw,
                             probs = seq(0.1, 0.9, length.out = K_counts),
                             names = FALSE)
prior_betas      <- prior_alpha / sort(pmax(mean_q, 10), decreasing = TRUE)
prior_T_counts   <- generate_random_T(K_counts, self = 0.85)

# ───── Gaussian priors (GC proportion) ───────────────────────────────────────
prior_sd         <- sd(gc_prop)
prior_means      <- sort(quantile(gc_prop,
                                  probs = seq(0.1, 0.9, length.out = K_gc),
                                  names = FALSE))
prior_T_gc       <- generate_random_T(K_gc, self = 0.85)

## ────────────────────────────────────────────────────────────────────────────
## 4.  TEST RUN of oHMMed (optional)
## ────────────────────────────────────────────────────────────────────────────
options(error = recover)    # drop into debugger on any unexpected error

set.seed(42)
fit_snv <- hmm_mcmc_gamma_poisson(
  data         = counts_raw,
  prior_T      = prior_T_counts,
  prior_betas  = prior_betas,
  prior_alpha  = prior_alpha,
  iter         = 8000,
  warmup       =   200,
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
  iter         = 8000,
  warmup       =   2000,
  thin         =    10,
  print_params = FALSE,
  verbose      = TRUE
)

## ────────────────────────────────────────────────────────────────────────────
## 5.  DIAGNOSTICS  +  SAVE EVERY FIGURE
## ────────────────────────────────────────────────────────────────────────────
output_dir <- "data/treated data/plots"          # <- where images will go
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## ───── 5a. built-in diagnostic plots ────────────────────────────────────────
# oHMMed::plot() returns a grid object; wrap it in a PNG device
png(file.path(output_dir, "diagnostics_snv_K4.png"),
    width = 1600, height = 1200, res = 150)
plot(fit_snv, main = "oHMMed diagnostics — SNV (K = 4)")
dev.off()

png(file.path(output_dir, "diagnostics_gc_K4.png"),
    width = 1600, height = 1200, res = 150)
plot(fit_gc,  main = "oHMMed diagnostics — GC (K = 4)")
dev.off()

## ───── 5b. ggmcmc trace plots ──────────────────────────────────────────────


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

## ────────────────────────────────────────────────────────────────────────────
## 7a.  SWEEP OVER K FOR GC PROPORTION  +  SAVE EVERY EXTRA PLOT
## ────────────────────────────────────────────────────────────────────────────
## ───── set K range for GC proportion ────────────────────────────────────────
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
    iter         = 8000, warmup = 2000, thin = 10,
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

## ───── log-likelihood summary figure ────────────────────────────────────────
loglik_gc_summary <- data.frame(
  K        = states,
  MeanLL   = vapply(fits_gc, function(f) mean(f$estimates$log_likelihood),
                    numeric(1)),
  MedianLL = vapply(fits_gc, function(f) median(f$estimates$log_likelihood),
                    numeric(1))
)


p_ll <- ggplot(loglik_gc_summary, aes(x = K)) +
          geom_line(aes(y = MeanLL, color = "Mean"),    linewidth = 1) +        # <- linewidth!
          geom_point(aes(y = MeanLL, color = "Mean"),   size = 2) +
          geom_line(aes(y = MedianLL, color = "Median"),  linetype = "dashed", linewidth = 1) +
          geom_point(aes(y = MedianLL, color = "Median"), size = 2) +
          scale_x_continuous(breaks = states) +
          scale_color_manual(values = c(Mean = "steelblue",   # teal
                                Median = "darkorange"))+ # orange
          labs(
            title = "Arabidopsis GC: mean / median log-likelihood vs. K",
            x     = "Number of hidden states (K)",
            y     = "Log-likelihood"
          ) +
          theme_bw()

ggsave(filename = file.path(output_dir, "gc_loglik_vs_K.png"),
       plot = p_ll, width = 9, height = 6, dpi = 150)


## ───── choose a K and extract the posterior states ───────────────────────────────
k5_idx  <- which(states == 5)                     # <- change this to choose K
fit_k5  <- fits_gc[[k5_idx]]
## get means
mu_draws <- fit_k5$estimates$mean   
mu_draws
## get sd 
fit_k5[["estimates"]][["sd"]]

win_state <- fit_k5$estimates$posterior_states     # length = n_windows

## counts and proportions
tbl_windows <- table(win_state)                    # counts per state
prop_windows <- prop.table(tbl_windows)            # proportions (0‒1)

prop_df <- data.frame(
  state      = as.integer(names(tbl_windows)),
  n_windows  = as.integer(tbl_windows),
  proportion = as.numeric(prop_windows)
)
prop_df

#### Genome-wide track  
gc_df$state <- fit_k5$estimates$posterior_states 
state_cols <- c("1"="#4C72B0",   # AT-rich
                "2"="#55A868",
                "3"="#E4BF44",
                "4"="#C35B72",
                "5"="#8172B2")   # GC-rich

bed <- gc_df |>
  dplyr::transmute(seqnames,
                   start,                  
                   end,                    
                   name      = paste0("S", state),
                   score     = 0,
                   strand    = "+",
                   thickStart = start,     
                   thickEnd   = end,
                   itemRgb    = state_cols[as.character(state)])   

genome_track <- ggplot(gc_df, aes(x = start/1e6,            
                  y = seqnames,                  # one row per chromosome
                  fill = factor(state))) +
  geom_tile(height = 0.8) +                 # 0.8 keeps thin horizontal gaps
  scale_fill_manual(values = state_cols,
                    name = "GC state") +
  scale_y_discrete(limits = paste0("Chr", 1:5)) +
  labs(x = "Chromosomal position (Mb)",
       y = NULL,
       title = "Five-state GC segmentation of the *A. thaliana* genome (100 kb windows)") +
  theme_bw() +
  theme(panel.spacing.y = unit(0.1, "lines"),
        axis.text.y     = element_text(size = 9),
        legend.position = "right")

ggsave(filename = file.path(output_dir, "Genome-wide track.png"),
       plot = genome_track, width = 9, height = 6, dpi = 150)

## ────────────────────────────────────────────────────────────────────────────
## 7b.  SWEEP OVER K  FOR SNV COUNTS   +   SAVE EVERY EXTRA PLOT
## ────────────────────────────────────────────────────────────────────────────
# ───── set K range for SNV counts ────────────────────────────────────────────
states_counts <- 2:9
fits_counts   <- vector("list", length(states_counts))

for (i in seq_along(states_counts)) {
  k <- states_counts[i]
  message("◆ Fitting SNV counts, K = ", k)
  set.seed(100 + k)
  
  # rebuild priors for this k
  mean_q_k        <- quantile(counts_raw,
                              probs = seq(0.1, 0.9, length.out = k),
                              names = FALSE)
  prior_betas_k   <- prior_alpha / sort(pmax(mean_q_k, 10), decreasing = TRUE)
  prior_T_counts_k <- generate_random_T(k, self = 0.85)
  
  
  # fit the gamma–Poisson HMM
  fits_counts[[i]] <- hmm_mcmc_gamma_poisson(
    data         = counts_raw,
    prior_T      = prior_T_counts_k,
    prior_betas  = prior_betas_k,
    prior_alpha  = prior_alpha,
    iter         = 8000,
    warmup       = 2000,
    thin         = 10,
    print_params = FALSE,
    verbose      = FALSE
  )
  
  # save diagnostic plot
  png(file.path(output_dir,
                sprintf("diagnostics_counts_K%02d.png", k)),
      width = 1600, height = 1200, res = 150)
  plot(fits_counts[[i]],
       main = sprintf("oHMMed diagnostics — SNV counts, K = %d", k))
  dev.off()
}

# ───── log-likelihood summary figure for SNV counts ──────────────────────
loglik_counts_summary <- data.frame(
  K        = states_counts,
  MeanLL   = vapply(fits_counts,
                    function(f) mean(f$estimates$log_likelihood),
                    numeric(1)),
  MedianLL = vapply(fits_counts,
                    function(f) median(f$estimates$log_likelihood),
                    numeric(1))
)

library(ggplot2)
p_ll_counts <- ggplot(loglik_counts_summary, aes(x = K)) +
  geom_line(aes(y = MeanLL,   color = "Mean"),    linewidth = 1) +
  geom_point(aes(y = MeanLL,   color = "Mean"),   size = 2) +
  geom_line(aes(y = MedianLL, color = "Median"),
            linetype = "dashed", linewidth = 1) +
  geom_point(aes(y = MedianLL, color = "Median"), size = 2) +
  scale_x_continuous(breaks = states_counts) +
  scale_color_manual(values = c(Mean = "steelblue",
                                Median = "darkorange")) +
  labs(
    title = "Breast-cancer SNV counts: mean / median log-likelihood vs. K",
    x     = "Number of hidden states (K)",
    y     = "Log-likelihood"
  ) +
  theme_bw()

ggsave(filename = file.path(output_dir, "counts_loglik_vs_K.png"),
       plot    = p_ll_counts,
       width   = 9,
       height  = 6,
       dpi     = 150)

# ───── choose a K and extract the posterior states ────────────────────
k <- 5 # <- change this to choose K
message("◆ Fitting SNV counts, K = ", k)
set.seed(100 + k)

# rebuild priors for this k
mean_q_k        <- quantile(counts_raw,
                            probs = seq(0.1, 0.9, length.out = k),
                            names = FALSE)
prior_betas_k   <- prior_alpha / sort(pmax(mean_q_k, 10), decreasing = TRUE)
prior_T_counts_k <- generate_random_T(k, self = 0.85)


# fit the gamma–Poisson HMM
fits_counts <- hmm_mcmc_gamma_poisson(
  data         = counts_raw,
  prior_T      = prior_T_counts_k,
  prior_betas  = prior_betas_k,
  prior_alpha  = prior_alpha,
  iter         = 8000,
  warmup       = 2000,
  thin         = 10,
  print_params = FALSE,
  verbose      = FALSE
)

snv_df$state <- fits_counts$estimates$posterior_states 

## choose a colour for every state  (low → high SNV burden)
state_cols_snv <- c("1" = "#D4B9DA",   # pale
                    "2" = "#C994C7",
                    "3" = "#DF65B0",
                    "4" = "#CE1256",
                    "5" = "#660033")   # dark

bed_snv <- snv_df |>
  dplyr::transmute(seqnames,
                   start,
                   end,
                   name       = paste0("S", state),
                   score      = 0,
                   strand     = "+",
                   thickStart = start,
                   thickEnd   = end,
                   itemRgb    = state_cols_snv[as.character(state)])

write.table(bed_snv,
            file = "data/treated data/breast_cancer/TCGA-BH-A201-01A-11D-A14K-09_oHMMed_K5_SNV.bed",
            sep  = "\t", quote = FALSE, row.names = FALSE)

## quick chromosomal heat-strip keeps the y-order in snv_df$seqnames
snv_track <- ggplot(snv_df,
                    aes(x   = start/1e6,     # Mb
                        y   = factor(seqnames, levels = rev(unique(seqnames))),
                        fill = factor(state))) +
  geom_tile(height = 0.8) +
  scale_fill_manual(values = state_cols_snv,
                    name   = "SNV state") +
  labs(x = "Chromosomal position (Mb)",
       y = NULL,
       title = "five-state SNV segmentation (100 kb windows)") +
  theme_bw() +
  theme(panel.spacing.y = unit(0.1, "lines"),
        axis.text.y     = element_text(size = 9),
        legend.position = "right")

ggsave(file.path(output_dir, "Genome-wide_SNV_track.png"),
       snv_track, width = 9, height = 6, dpi = 150)