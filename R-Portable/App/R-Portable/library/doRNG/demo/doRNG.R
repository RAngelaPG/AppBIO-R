library(doRNG)

## ---- doParallel example ------------------------------------------------
if (requireNamespace("doParallel", quietly = TRUE)) {
  
  # PSOCK/SNOW-like backend (works on all OS)
  cl <- NULL
  on.exit({
    foreach::registerDoSEQ()
    if (!is.null(cl)) parallel::stopCluster(cl)
  }, add = TRUE)
  
  cl <- parallel::makeCluster(2)
  doParallel::registerDoParallel(cl)
  
  r1 <- foreach::foreach(i = 1:4, .options.RNG = 1234) %dorng% { runif(1) }
  r2 <- foreach::foreach(i = 1:4, .options.RNG = 1234) %dorng% { runif(1) }
  identical(r1, r2)
  attr(r1, "rng")
  
  set.seed(1234)
  s1 <- foreach::foreach(i = 1:4) %dorng% { runif(1) }
  s2 <- foreach::foreach(i = 1:4) %dorng% { runif(1) }
  identical(s1, s2)
  
  set.seed(1234)
  s1.2 <- foreach::foreach(i = 1:4) %dorng% { runif(1) }
  s2.2 <- foreach::foreach(i = 1:4) %dorng% { runif(1) }
  identical(s1, s1.2) && identical(s2, s2.2)
  identical(r1, s1)
  
  # Optional fork backend (Unix only) - keep it separate
  if (.Platform$OS.type == "unix") {
    doParallel::registerDoParallel(cores = 2)
    
    r1 <- foreach::foreach(i = 1:4, .options.RNG = 1234) %dorng% { runif(1) }
    r2 <- foreach::foreach(i = 1:4, .options.RNG = 1234) %dorng% { runif(1) }
    identical(r1, r2)
    
    foreach::registerDoSEQ()
  }
  
} else {
  message("Skipping doParallel demo: package 'doParallel' is not installed.")
}

## ---- doMPI example -----------------------------------------------------
if (requireNamespace("doMPI", quietly = TRUE)) {
  
  cl_mpi <- NULL
  on.exit({
    foreach::registerDoSEQ()
    if (!is.null(cl_mpi)) doMPI::closeCluster(cl_mpi)
  }, add = TRUE)
  
  cl_mpi <- doMPI::startMPIcluster(2)
  doMPI::registerDoMPI(cl_mpi)
  
  s1 <- foreach::foreach(i = 1:4, .options.RNG = 1234) %dorng% { runif(1) }
  s2 <- foreach::foreach(i = 1:4, .options.RNG = 1234) %dorng% { runif(1) }
  identical(s1, s2)
  
} else {
  message("Skipping doMPI demo: package 'doMPI' is not installed.")
}
