## ----style, echo=FALSE, results='hide', message=FALSE-------------------------
library(BiocStyle)
library(knitr)
opts_chunk$set(error = FALSE, message = FALSE, warning = FALSE)
opts_chunk$set(fig.asp = 1)

args <- commandArgs(trailingOnly = TRUE)
print(args)
transformpath = args[1]
predcluster = args[2]
n_clusters = as.integer(args[3])

suppressPackageStartupMessages(library(Melissa)) # Load package

## ----load_data
source("data_process.R")
cat("Number of cells: ", length(dt_obj$met))
cat("Number of genomic regions: ", length(dt_obj$met[[1]]) )
library(BPRMeth)
basis_obj <- create_rbf_object(M = 4)
basis_obj
set.seed(15)
# Partition to training and test set
dt_obj <- partition_dataset(dt_obj = dt_obj, data_train_prcg = 0.2,
                            region_train_prcg = 1, cpg_train_prcg = 0.4, 
                            is_synth = TRUE)

## ----run_melissa--------------------------------------------------------------
set.seed(15)
melissa_obj <- melissa(X = dt_obj$met, K = n_clusters, basis = basis_obj,
                       vb_max_iter = 20, vb_init_nstart = 1, 
                       is_parallel = FALSE)

## ----summary_mixing_proportions-----------------------------------------------
melissa_obj$pi_k

## ----summary_responsibilities-------------------------------------------------
print(melissa_obj$r_nk)

write.table(melissa_obj$r_nk, file = predcluster, sep = '\t', row.names = TRUE, col.names = TRUE)
sessionInfo()
