## ----missing-response, message=FALSE------------------------------------------
library(breedR)

N <- 1e3
x <- rep(1:4, each = N/4)
dat <- data.frame(y = x + rnorm(N),
                  x = factor(letters[x]))
dat$y[1] <- NA
head(dat)
res <- remlf90(y ~ x, data = dat)

## The predicted phenotype for y[1] is the estimated effect
## of the corresponding level of x
fitted(res)[1] == fixef(res)$x['a']

## ----missing-fixed, message = FALSE, error = TRUE-----------------------------
N <- 1e3
x <- rep(1:4, each = N/4)
dat <- data.frame(y = x + rnorm(N),
                  x = factor(letters[x]))
dat$x[c(1, 3, 5)] <- NA
head(dat)
res <- remlf90(y ~ x, data = dat)

## ----missing-fixed-regression, message = FALSE, error = TRUE------------------
N <- 1e3
x <- runif(N)
dat <- data.frame(y = 1 + 2*x + rnorm(N),
                  x = x)
dat$x[c(1, 3, 5)] <- NA
head(dat)
res <- remlf90(y ~ x, data = dat)

## ----missing-diagonal, message = FALSE----------------------------------------
N <- 1e3
N.blk <- 20
blk.effects <- rnorm(N.blk, sd = 2)
blk.idx <- sample(seq_len(N.blk), N, replace = TRUE)
dat <- data.frame(y = 1 + blk.effects[blk.idx] + rnorm(N),
                  blk = factor(blk.idx))
dat$blk[1] <- NA
head(dat)

res <- remlf90(y ~ 1, random = ~ blk, data = dat)

sum(model.matrix(res)$blk[1,])

## ----missing-diagonal-residual------------------------------------------------
fitted(res)[1] == fixef(res)$Intercept[1]

## ----missing-block, message = FALSE-------------------------------------------
coord <- expand.grid(row = 1:20, col = 1:50)
res <- remlf90(y ~ 1,
               spatial = list(model = 'blocks',
                              coord = coord,
                              id    = 'blk'), 
               data = dat)

c(sum(model.matrix(res)$spatial[1,]) == 0,
fitted(res)[1] == fixef(res)$Intercept[1])

## ----sample-residual-function, echo = FALSE-----------------------------------
sample_first_residual <- function(N = 1e3, N.blk = 20) {
  blk.effects <- rnorm(N.blk, sd = 2)
  blk.idx <- sample(seq_len(N.blk), N, replace = TRUE)
  dat <- data.frame(y = 1 + blk.effects[blk.idx] + rnorm(N),
                    blk = factor(blk.idx))
  dat$blk[1] <- NA
  res <- suppressMessages(remlf90(y ~ 1, random = ~ blk, data = dat))
  return(unname(resid(res)[1]))
}

## ----variance-missing-residuals-----------------------------------------------
resid_sample <- replicate(1e2, sample_first_residual())
var(resid_sample)


## ----missing-add_animal, message = FALSE--------------------------------------
dat <- breedR.sample.phenotype(
  fixed = c(mu = 10, x = 2),
  genetic = list(model    = 'add_animal',
                 Nparents = c(10, 10),
                 sigma2_a = 2,
                 check.factorial = FALSE),
  N = 1e3)
head(dat)

res <- remlf90(phenotype ~ 1 + X.x,
               genetic = list(model = 'add_animal',
                              pedigree = dat[, 1:3],
                              id    = 'self'),
               data = dat)

str(ranef(res)$genetic)

## ----missing-coordinates, message = FALSE-------------------------------------
dat <- breedR.sample.phenotype(
  fixed = c(mu = 10, x = 2),
  spatial = list(model     = 'AR',
                 grid.size = c(10, 5),
                 rho       = c(.2, .8),
                 sigma2_s  = 1)
)
dat$Var1[1] <- NA
head(dat)

res <- remlf90(phenotype ~ 1 + X.x,
               spatial = list(model = 'AR',
                              coord = dat[, c('Var1', 'Var2')],
                              rho   = c(0.2, 0.8)),
               data = dat)

sum(model.matrix(res)$spatial[1,])

