library(mgcv)
data("globulus")

object <- tb_spec <- s(self, dad, mum, bs = "pd")
data <- globulus
knots <- NULL

smooth.construct.pd.smooth.spec <- function(object,data,knots) {

  ped <- build_pedigree(1:3, data = data[object$term])
  
  add_animal <- additive_genetic_animal(ped, data[[object$term[1]]])

  object$X <- as.matrix(add_animal$incidence.matrix)
  object$S <- list(as.matrix(add_animal$structure.matrix))
  object$rank <- ncol(object$S)  # full rank
  object$null.space.dim <- 0
    
  class(object) <- "pedigree.smooth"
  object
}

res.animal <- remlf90(
  fixed  = phe_X ~ gg,
  genetic = list(
    model = 'add_animal', 
    pedigree = globulus[, 1:3],
    id = 'self'), 
  data = globulus)

resmgcv <- gam(phe_X ~ gg + s(self, dad, mum, bs = "pd"), data = globulus)

# Error in gam(phe_X ~ gg + s(self, dad, mum, bs = "pd"), data = globulus) : 
#   Model has more coefficients than data

## Apparently, gam() checks that the model matrix is longer than it 
## is wide. But this is not the case for a additive-genetic effect.
## 

# Predict.matrix.tr.smooth<-function(object,data)
#   ## prediction method function for the `tr' smooth class
# { x <- data[[object$term]]
# x <- x - object$x.shift # stabilizing shift
# m <- object$m;     # spline order (3=cubic)
# k<-object$knots    # knot locations
# nk<-length(k)      # number of knots
# X<-matrix(0,length(x),object$bs.dim)
# for (i in 1:(m+1)) X[,i] <- x^(i-1)
# for (i in 1:nk) X[,i+m+1] <- (x-k[i])^m*as.numeric(x>k[i])
# X # return the prediction matrix
# }
# 
# # an example, using the new class....
# require(mgcv)
# set.seed(100)
# dat <- gamSim(1,n=400,scale=2)
# b<-gam(y~s(x0,bs="tr",m=2)+s(x1,bs="ps",m=c(1,3))+
#          s(x2,bs="tr",m=3)+s(x3,bs="tr",m=2),data=dat)
# plot(b,pages=1)
# b<-gamm(y~s(x0,bs="tr",m=2)+s(x1,bs="ps",m=c(1,3))+
#           s(x2,bs="tr",m=3)+s(x3,bs="tr",m=2),data=dat)
# plot(b$gam,pages=1)
# # another example using tensor products of the new class
# dat <- gamSim(2,n=400,scale=.1)$data
# b <- gam(y~te(x,z,bs=c("tr","tr"),m=c(2,2)),data=dat)
# vis.gam(b)
