library(oHMMed)

df <- read.table("~/Desktop/bioinf/M2/methods bionfo /GenomeAnnotations/HumanAnnotation.txt",
                header = TRUE, 
                sep="")


gc_values <- df$GC_PROPORTION

K <- 5 
#prior_means <- seq(0.2869171,0.5474,length.out = K)
prior_means <- seq(min(gc_values), max(gc_values), length.out = 5)
prior_sd <-  2.5 / 5
prior_T <- generate_random_T(K)
iter <- 1500

model <- hmm_mcmc_normal(
  data = gc_values,
  prior_T = prior_T,
  prior_means = prior_means,
  prior_sd = prior_sd,
  iter = iter,
  warmup = floor(iter * 0.2),
  print_params = FALSE,
  verbose = TRUE
)


summary(model)
plot(model)


