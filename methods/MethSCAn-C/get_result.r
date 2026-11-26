library(tidyverse)
library(irlba)

args <- commandArgs(trailingOnly=TRUE)
output_VMR_matrix = args[1]

prcomp_iterative <- function(x, n=10, n_iter=50, min_gain=0.001, ...) {
  mse <- rep(NA, n_iter)
  na_loc <- is.na(x)
  x[na_loc] = 0  # zero is our first guess

  for (i in 1:n_iter) {
    prev_imp <- x[na_loc]  # what we imputed in the previous round
    # PCA on the imputed matrix
    pr <- prcomp_irlba(x, center = F, scale. = F, n = n, ...)
    # impute missing values with PCA
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    # compare our new imputed values to the ones from the previous round
    mse[i] = mean((prev_imp - new_imp) ^ 2)
    # if the values didn't change a lot, terminate the iteration
    gain <- mse[i] / max(mse, na.rm = T)
    if (gain < min_gain) {
      message(paste(c("\n\nTerminated after ", i, " iterations.")))
      break
    }
  }
  pr$mse_iter <- mse[1:i]
  pr
}

meth_mtx <- read.csv(file.path(output_VMR_matrix, "mean_shrunken_residuals.csv.gz"), row.names=1) %>%
  as.matrix()

pca <- meth_mtx %>%
  scale(center = T, scale = F) %>%
  prcomp_iterative(n = 15)  # increase this value to e.g. 15 for real data sets

pca_tbl <- as_tibble(pca$x) %>% 
  add_column(cell=rownames(meth_mtx))

pca_loading <- as_tibble(pca$rotation)

readr::write_csv(
    pca_loading, 
    file = file.path(output_VMR_matrix, "pca_loading_mtx.csv")
)

readr::write_csv(
    pca_tbl, 
    file = file.path(output_VMR_matrix, "pca_results_mtx.csv")
)
