############################################################
# Clean environment and set working directory
############################################################

rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


############################################################
# Packages
############################################################

library(Rcpp)
library(RcppEigen)
library(stringr)
library(parallel)
library(pbapply)


############################################################
# Source files
############################################################

source("Gen_determinData.R")
source("test.R")


############################################################
# Simulation settings
############################################################

n = c(400, 800, 1200)
n = c(1200) ####
MCMC = 500

alpha = c(0.01, 0.05, 0.1, 0.2)


############################################################
# Create output directory if necessary
############################################################

if (!dir.exists("result")) {
  dir.create("result")
}


############################################################
# Parallel settings
############################################################

# 建议先用4核
# 如果内存足够，可以改成6或者8
n_cores = 4

cat("Detected CPU cores:", detectCores(), "\n")
cat("Using", n_cores, "cores for parallel computation.\n\n")

cl = makeCluster(n_cores)


############################################################
# Set working directory for every worker
############################################################

work_dir = getwd()

clusterExport(
  cl,
  varlist = "work_dir",
  envir = .GlobalEnv
)


############################################################
# Load functions/packages on every worker
############################################################

clusterEvalQ(cl, {
  
  setwd(work_dir)
  
  # test.R 会进一步 source Gen_randomdata.R 和 Wu.R
  source("Gen_determinData.R")
  source("test.R")
  
  NULL
})


############################################################
# Reproducible parallel random-number streams
############################################################

clusterSetRNGStream(cl, iseed = 123)


############################################################
# Progress bar
############################################################

# 显示百分比、运行时间和预计剩余时间
pboptions(type = "timer")


############################################################
# Function for one sample size n
############################################################

f = function(n) {
  
  N = 100
  type = "AR"
  
  cat("\n")
  cat("====================================================\n")
  cat("Starting simulation\n")
  cat("T =", n, "\n")
  cat("N =", N, "\n")
  cat("MCMC =", MCMC, "\n")
  cat("Number of cores =", n_cores, "\n")
  cat("====================================================\n\n")
  
  start_time = Sys.time()
  
  
  ##########################################################
  # Deterministic part of the data-generating mechanism
  ##########################################################
  
  Data = Gen_determinData(N)
  
  
  ##########################################################
  # Parallel Monte Carlo
  ##########################################################
  
  raw_list = pblapply(
    X = 1:MCMC,
    FUN = test,
    
    N = N,
    n = n,
    alpha = alpha,
    Data = Data,
    type = type,
    
    cl = cl
  )
  
  
  ##########################################################
  # Convert list to matrix
  #
  # Each column corresponds to one Monte Carlo replication
  ##########################################################
  
  raw_result = do.call(cbind, raw_list)
  
  
  ##########################################################
  # Average across Monte Carlo replications
  ##########################################################
  
  avg_result = rowMeans(
    raw_result,
    na.rm = TRUE
  )
  
  
  ##########################################################
  # Arrange results
  #
  # 6 scenarios:
  # Homo Normal / Uniform / Laplace
  # Hetero Normal / Uniform / Laplace
  #
  # Each scenario contains:
  # 4 full-domain coverage probabilities
  # 4 full-domain volumes
  # 4 interior coverage probabilities
  # 4 interior volumes
  ##########################################################
  
  result = matrix(
    avg_result,
    nrow = 6,
    ncol = 4 * length(alpha),
    byrow = TRUE
  )
  
  
  ##########################################################
  # Row names
  ##########################################################
  
  rownames(result) = c(
    "Homo_Normal",
    "Homo_Uniform",
    "Homo_Laplace",
    "Hetero_Normal",
    "Hetero_Uniform",
    "Hetero_Laplace"
  )
  
  
  ##########################################################
  # Column names
  ##########################################################
  
  colnames(result) = c(
    
    paste0(
      "Coverage_full_alpha_",
      alpha
    ),
    
    paste0(
      "Volume_full_alpha_",
      alpha
    ),
    
    paste0(
      "Coverage_interior_alpha_",
      alpha
    ),
    
    paste0(
      "Volume_interior_alpha_",
      alpha
    )
  )
  
  
  ##########################################################
  # Print result
  ##########################################################
  
  cat("\n")
  cat("Finished T =", n, "\n\n")
  
  print(result)
  
  
  ##########################################################
  # Save result
  ##########################################################
  
  output_file = str_c(
    "result/Wu_T=",
    n,
    "N=",
    N,
    ".csv"
  )
  
  write.csv(
    result,
    file = output_file
  )
  
  cat("\nSaved to:", output_file, "\n")
  
  
  ##########################################################
  # Running time
  ##########################################################
  
  end_time = Sys.time()
  
  cat(
    "Running time for T =",
    n,
    ":",
    as.numeric(
      difftime(
        end_time,
        start_time,
        units = "mins"
      )
    ),
    "minutes\n"
  )
  
  cat("====================================================\n\n")
  
  
  return(result)
}


############################################################
# Run T = 400, 800, 1200
############################################################

overall_start = Sys.time()


A = lapply(
  n,
  f
)


############################################################
# Stop parallel workers
############################################################

stopCluster(cl)


############################################################
# Total running time
############################################################

overall_end = Sys.time()

cat("\n")
cat("====================================================\n")
cat("All simulations finished.\n")

cat(
  "Total running time:",
  as.numeric(
    difftime(
      overall_end,
      overall_start,
      units = "hours"
    )
  ),
  "hours\n"
)

cat("====================================================\n")