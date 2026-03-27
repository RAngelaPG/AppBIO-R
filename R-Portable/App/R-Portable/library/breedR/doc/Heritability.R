## ----setup, include = FALSE---------------------------------------------------
library(breedR)
library(ggplot2)

## ----fit-additive-genetic-ai, message = FALSE---------------------------------
res.animal <- remlf90(fixed  = phe_X ~ 1,
                      random = ~ gg,
                      genetic = list(model = 'add_animal', 
                                     pedigree = globulus[, 1:3],
                                     id = 'self'), 
                      data = globulus)
summary(res.animal)

## ----globulus-se_covar_function, message = FALSE------------------------------
globulus <- transform(globulus,
                      fam = factor(mum, exclude = 0))

h2fml <- '4*G_3_3_1_1/(G_3_3_1_1+R_1_1)'

res.fml <- remlf90(fixed  = phe_X ~ gg,
                   random = ~ bl + fam,
                   progsf90.options = paste('se_covar_function h2', h2fml),
                   data = globulus)
summary(res.fml)

## ----extract-invAI------------------------------------------------------------
(invAI <- res.fml$reml$invAI)

## ----sqrt-invAI---------------------------------------------------------------
sqrt(diag(invAI))

## ----deltamethod-msm----------------------------------------------------------
if (require(msm, quietly = TRUE)) {
  
  h2hat <- with(as.list(res.fml$var[, "Estimated variances"]),
                4*fam/(fam+Residual))
  h2se <- deltamethod(~ 4*x2/(x2+x3),
                      res.fml$var[, "Estimated variances"],
                      invAI)
  cat(paste0('Heritability (delta method): ', 
             round(h2hat, 2), ' (', 
             round(h2se, 2), ')'))
}

## ----simple-model-fit, message = FALSE----------------------------------------
globulus <- transform(globulus,
                      fam = factor(mum, exclude = 0))

res <- remlf90(fixed  = phe_X ~ 1,
               random = ~ bl + fam,
               data = globulus)

## ----resampling-function------------------------------------------------------

resample_globulus <- function(N = 1, fit = res) {
  
  dat <- breedR.sample.phenotype(
    fixed = 1,
    random = list(bl = list(nlevels = 15,
                            sigma2 = fit$var['bl', 1]),
                  fam = list(nlevels = 63,
                             sigma2 = fit$var['fam', 1])),
    residual.variance = fit$var['Residual', 1],
    N = nrow(model.frame(fit))
  )
  
  # rename variables to keep the original names
  names(dat) <- gsub("rnd\\.", "", names(dat))
  
  return(dat)
}

## ----simulate-variances, message = FALSE--------------------------------------

sim_target <- function(dat = resample_globulus()) {

  ## Fit the original model
  ## avoiding messages about initial variances spec
  res <- suppressMessages(
    remlf90(fixed  = phenotype ~ 1,
            random = ~ bl + fam,
            data = dat)
  )
  
  ## Extract point estimates of all variances
  ## and point estimate of heritability
  ans <- res$var[, 'Estimated variances']
  ans <- c(ans, h2 = unname(4*ans['fam']/sum(ans)))

  return(ans)
}


## ----simulator-example--------------------------------------------------------
sim_target()
sim_target()
sim_target()

## ----sampling-distribution-variances, fig.cap = 'Approximate sampling distribution of variance components', message = FALSE----
empirical_dist <- as.data.frame(t(replicate(500, sim_target())))

if (require(tidyr, quietly = TRUE)) {

  plotdat <- gather(empirical_dist[, 1:3], 'variable', 'value')
  qplot(value, data = plotdat, facets = ~ variable)
}


## ----sampling-distribution-heritability, fig.cap = 'Approximate sampling distribution of heritability', message = FALSE----

qplot(h2, data = empirical_dist, xlim = c(0, 1))
  

## ----standard-errors, results='hold'------------------------------------------
cat('Standard Errors:\n')
round(sapply(empirical_dist, sd), 2)

## ----confidence-intervals, results='hold'-------------------------------------
cat('95% Confidence Intervals:\n')
round(sapply(empirical_dist, quantile, probs = c(.025, .975)), 2)

