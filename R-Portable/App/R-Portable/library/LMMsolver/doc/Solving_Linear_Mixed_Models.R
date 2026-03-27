## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dim = c(7, 4)
)

## ----setup--------------------------------------------------------------------
library(LMMsolver)
library(ggplot2)

## ----linearFun----------------------------------------------------------------
  f1 <- function(x) { 0.6 + 0.7*x}

## ----linearFunSim-------------------------------------------------------------
set.seed(2016)
n <- 25
x <- seq(0, 1, length = n)
sigma2e <- 0.04
y <- f1(x) + rnorm(n, sd = sqrt(sigma2e))
dat1 <- data.frame(x = x, y = y)

## ----fit1---------------------------------------------------------------------
obj1 <- LMMsolve(fixed = y ~ x, data = dat1)

## ----predictlin---------------------------------------------------------------
newdat <- data.frame(x = seq(0, 1, length = 300))
pred1 <- predict(obj1, newdata = newdat, se.fit = TRUE)
# adding the true values for comparison
pred1$y_true <- f1(pred1$x)

## ----ggplotLinear-------------------------------------------------------------
ggplot(data = dat1, aes(x = x, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred1, aes(y=y_true), color = "red", 
            linewidth = 1, linetype = "dashed") +
  geom_line(data = pred1, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data = pred1, aes(x=x,ymin = ypred-2*se, ymax = ypred+2*se),
              alpha = 0.2, inherit.aes = FALSE) + 
  theme_bw()

## ----nonlinearFun-------------------------------------------------------------
g <- function(x) { 0.3 + 0.4*x + 0.2*sin(20*x) }

## ----simDat2------------------------------------------------------------------
set.seed(12)
n <- 150
x <- seq(0, 1, length = n)
sigma2e <- 0.04
y <- g(x) + rnorm(n, sd = sqrt(sigma2e))
dat2 <- data.frame(x, y)

## ----nonlinearFit-------------------------------------------------------------
obj2 <- LMMsolve(fixed = y ~ 1, 
                 spline = ~spl1D(x, nseg = 50), 
                 data = dat2)

## ----summary_nonlinear--------------------------------------------------------
summary(obj2)

## ----predict_non--------------------------------------------------------------
newdat <- data.frame(x = seq(0, 1, length = 300))
pred2 <- predict(obj2, newdata = newdat, se.fit = TRUE)
pred2$y_true <- g(pred2$x)

ggplot(data = dat2, aes(x = x, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred2, aes(y = y_true), color = "red", 
            linewidth = 1, linetype ="dashed") +
  geom_line(data = pred2, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred2, aes(x=x, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  theme_bw() 

## ----predict_non_derivative---------------------------------------------------
dg_dx <- function(x) { 0.4 + 4*cos(20*x) }
pred2_dx <- predict(obj2, newdata = newdat, se.fit = TRUE, deriv="x")
pred2_dx$y_true <- dg_dx(pred2_dx$x)

ggplot(data = pred2_dx) +
  geom_line(aes(x=x, y = y_true), color = "red", 
            linewidth = 1, linetype ="dashed") +
  geom_line(aes(x = x, y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(aes(x=x, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  geom_hline(yintercept=0.0, linetype = "dashed") + ylab('dy/dx') + 
  theme_bw() 

## ----twoExperiments-----------------------------------------------------------
set.seed(1234)
nA <-  50
nB <- 100
mu_A <-  0.10
mu_B <- -0.10
sigma2e_A <- 0.04
sigma2e_B <- 0.10

x1 <- runif(n = nA)
x2 <- runif(n = nB)
y1 <- g(x1) + rnorm(nA, sd = sqrt(sigma2e_A)) + mu_A
y2 <- g(x2) + rnorm(nB, sd = sqrt(sigma2e_B)) + mu_B
Experiment <- as.factor(c(rep("A", nA), rep("B", nB)))
dat4 <- data.frame(x = c(x1, x2), y = c(y1,y2), Experiment = Experiment)

## ----boxplot------------------------------------------------------------------
ggplot(dat4, aes(x = Experiment, y = y, fill = Experiment)) +  
  geom_boxplot() + 
  geom_point(position = position_jitterdodge(), alpha = 0.3) 

## ----twoExpSolve--------------------------------------------------------------
obj4 <- LMMsolve(fixed= y ~ 1, 
                 spline = ~spl1D(x, nseg = 50, xlim = c(0,1)),
                 random = ~Experiment,
                 residual = ~Experiment,
                 data = dat4)

## ----summaryExp---------------------------------------------------------------
summary(obj4)

## ----twoExpPredict------------------------------------------------------------
newdat <- data.frame(x=seq(0, 1, length = 300))
pred4 <- predict(obj4, newdata = newdat, se.fit = TRUE)
pred4$y_true <- g(pred4$x)
ggplot(data = dat4, aes(x = x, y = y, colour = Experiment)) +
  geom_point(size = 1.5) +
  geom_line(data = pred4, aes(y = y_true), color="red", 
            linewidth = 1, linetype = "dashed") +
  geom_line(data = pred4, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data = pred4, aes(x = x,ymin = ypred-2*se, ymax = ypred+2*se),
              alpha = 0.2, inherit.aes = FALSE) + 
  theme_bw()

## ----coefExp------------------------------------------------------------------
coef(obj4)$Experiment

## ----USprecip data------------------------------------------------------------
## Get precipitation data from spam
data(USprecip, package = "spam")

## Only use observed data
USprecip <- as.data.frame(USprecip)
USprecip <- USprecip[USprecip$infill == 1, ]

## ----runobj3------------------------------------------------------------------
obj5 <- LMMsolve(fixed = anomaly ~ 1,
                 spline = ~spl2D(x1 = lon, x2 = lat, nseg = c(41, 41)),
                 data = USprecip)

## ----ED_USprecip--------------------------------------------------------------
summary(obj5)

## ----pred_USprecip------------------------------------------------------------
lon_range <- range(USprecip$lon)
lat_range <- range(USprecip$lat)
newdat <- expand.grid(lon = seq(lon_range[1], lon_range[2], length = 200),
                      lat = seq(lat_range[1], lat_range[2], length = 300))
plotDat5 <- predict(obj5, newdata = newdat)

## ----Plot_USprecip------------------------------------------------------------
plotDat5 <- sf::st_as_sf(plotDat5, coords = c("lon", "lat"))
usa <- sf::st_as_sf(maps::map("usa", regions = "main", plot = FALSE))
sf::st_crs(usa) <- sf::st_crs(plotDat5)
intersection <- sf::st_intersects(plotDat5, usa)
plotDat5 <- plotDat5[!is.na(as.numeric(intersection)), ]

ggplot(usa) + 
  geom_sf(color = NA) +
  geom_tile(data = plotDat5, 
            mapping = aes(geometry = geometry, fill = ypred), 
            linewidth = 0,
            stat = "sf_coordinates") +
  scale_fill_gradientn(colors = topo.colors(100))+
  labs(title = "Precipitation (anomaly)", 
       x = "Longitude", y = "Latitude") +
  coord_sf() +
  theme(panel.grid = element_blank())

## ----SeaSurfaceTemp-----------------------------------------------------------
data(SeaSurfaceTemp)
head(SeaSurfaceTemp, 3)
table(SeaSurfaceTemp$type)

## ----convert_and_split--------------------------------------------------------
# convert from Kelvin to Celsius
df <- SeaSurfaceTemp
df$sst <- df$sst - 273.15
### split in training and test set
df_train <- df[df$type == "train", ]
df_test <- df[df$type == "test", ]


## ----SST_raw_data-------------------------------------------------------------
nasa_palette <- c(
  "#03006d","#02008f","#0000b6","#0001ef","#0000f6","#0428f6","#0b53f7",
  "#0f81f3","#18b1f5","#1ff0f7","#27fada","#3efaa3","#5dfc7b","#85fd4e",
  "#aefc2a","#e9fc0d","#f6da0c","#f5a009","#f6780a","#f34a09","#f2210a",
  "#f50008","#d90009","#a80109","#730005"
)

map_layer <- geom_map(
  data = map_data("world"), map = map_data("world"),
  aes(group = group, map_id = region),
  fill = "black", colour = "white", linewidth = 0.1
)

# Brazil-Malvinas confluence zone
BM_box <- cbind(lon = c(-60, -48), lat = c(-50, -35))

ggplot() +
  scale_colour_gradientn(colours = nasa_palette, name = expression(degree*C)) +
  xlab("Longitude (deg)") + ylab("Latitude (deg)") +
  map_layer + xlim(BM_box[, "lon"]) + ylim(BM_box[, "lat"]) + theme_bw() +
  coord_fixed(expand = FALSE) +
  geom_point(data = df_train, aes(lon, lat, colour = sst), size=0.5)

## ----SeaTempAnalysis----------------------------------------------------------
obj6 <- LMMsolve(fixed = sst ~ 1, 
                 spline = ~spl2D(lon, lat, nseg = c(70, 70),
                                 x1lim = BM_box[, "lon"], x2lim = BM_box[, "lat"]),
                 data = df_train, tolerance = 1.0e-1)
summary(obj6)

## ----predictions_grid---------------------------------------------------------
lon_range <- BM_box[, "lon"]
lat_range <- BM_box[, "lat"]
newdat <- expand.grid(lon = seq(lon_range[1], lon_range[2], length = 200),
                      lat = seq(lat_range[1], lat_range[2], length = 200))

pred_grid <- predict(obj6, newdata = newdat, se.fit=TRUE)
pred_grid <- pred_grid[pred_grid$se<5, ]

## Plot predictions on a grid
ggplot(pred_grid) +
  geom_tile(aes(x = lon, y = lat, fill = ypred)) +
  scale_fill_gradientn(colours = nasa_palette) +
  labs(
    fill = "pred.",
    x = "Longitude (deg)", y = "Latitude (deg)"
  ) +
  map_layer +
  theme_bw() +
  coord_fixed(expand = FALSE, xlim = BM_box[, "lon"], ylim = BM_box[, "lat"]) +
  scale_x_continuous(breaks = c(-58, -54, -50))

## ----seSST--------------------------------------------------------------------
 ## Plot standard error
ggplot(pred_grid) +
   geom_raster(aes(x = lon, y = lat, fill = se)) +
   scale_fill_distiller(palette = "BrBG", direction = -1) +
   labs( fill = "s.e.", x = "Longitude (deg)", y = "Latitude (deg)") +
   map_layer +
   theme_bw() +
   coord_fixed(expand = FALSE, xlim = c(-60, -48), ylim = c(-50, -35)) +
   scale_x_continuous(breaks = c(-58, -54, -50))

## ----test_set-----------------------------------------------------------------
pred_test <- predict(obj6, newdata = df_test)
ggplot(pred_test, aes(x = sst,y = ypred)) + geom_point() +
     xlab("observed SST (Celsius)") + ylab("predicted SST (Celsius)") +
     geom_abline(intercept=0,slope=1,col='red') + theme_bw()

## ----RMSE_test----------------------------------------------------------------
Y <- (pred_test$sst - pred_test$ypred)^2
RMSE <- sqrt(mean(Y))
round(RMSE, 2)

## ----poissonSim---------------------------------------------------------------
set.seed(1234)
n <- 150
x <- seq(0, 1, length=n)
fun_lambda <- function(x) { 4 + 3*x + 4*sin(7*x) }
x <- seq(0, 1, length = n)
y <- rpois(n = n, lambda = fun_lambda(x)) 
dat3 <- data.frame(x = x, y = y)

## ----Poisson_LMMsolver--------------------------------------------------------
obj3 <- LMMsolve(fixed = y ~ 1, 
                spline = ~spl1D(x, nseg = 50), 
                family = poisson(), 
                data = dat3)
summary(obj3)

## ----predict_poisson----------------------------------------------------------
newdat <- data.frame(x = seq(0, 1, length = 300))
pred3 <- predict(obj3, newdata = newdat, se.fit = TRUE)
pred3$y_true <- fun_lambda(pred3$x)

ggplot(data = dat3, aes(x = x, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred3, aes(y = y_true), color = "red", 
            linewidth = 1, linetype ="dashed") +
  geom_line(data = pred3, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred3, aes(x=x, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  theme_bw() 

## ----binomial_sim-------------------------------------------------------------
set.seed(1234)
n <- 100
sz <- 10

fun_prob <- function(x) { 0.5 + 0.4*sin(2*pi*x) }

x <- seq(0, 1, length=n)
nsucces <- sapply(x, FUN=function(x) {
                  rbinom(n=1, size = sz, prob = fun_prob(x))
                })
dat <- data.frame(x = x, succes = nsucces, 
                         failure= sz - nsucces)
head(dat, 3)

## ----binomial_solve-----------------------------------------------------------
obj3 <- LMMsolve(fixed = cbind(succes, failure) ~ 1,
                 spline = ~spl1D(x, nseg = 50),
                 family = binomial(),
                 data = dat)
summary(obj3)

## ----binomial_predict---------------------------------------------------------
newdat <- data.frame(x = seq(0, 1, by=0.01))
pred3 <- predict(obj3, newdata = newdat, se.fit=TRUE)

## ----binomial_plot------------------------------------------------------------
pred3$y_true <- fun_prob(pred3$x)
dat$y <- dat$succes/sz

ggplot(data = dat, aes(x = x, y = y)) +
  geom_point(col = "black", size = 1.5) +
  geom_line(data = pred3, aes(y = y_true), color = "red",
            linewidth = 1, linetype = "dashed") +
  geom_line(data = pred3, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data= pred3, aes(x=x, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +
  theme_bw()

## ----multinomialmodel---------------------------------------------------------
k <- 4
mu <- c(0.1, 0.4, 0.6, 0.9)
names(mu) <- LETTERS[1:k]

nonlinear <- function(x, mu) {
  z <- sapply(mu, function(mu) { exp(-8*sin(pi*(x-mu))^2)})
  # normalize to sum equal to one
  z <- z/sum(z)
  return(z)
}

## ----simMultinomial-----------------------------------------------------------
x <- runif(n, 0, 1)   
sz <- 10 
multiNom <- t(sapply(x, FUN=
                      function(x) {
                        rmultinom(n=1, size=sz, prob = nonlinear(x,mu))
                      } ))
colnames(multiNom) <- names(mu)
dat <- data.frame(x, multiNom)
head(dat, 3)

## ----multinomialfit-----------------------------------------------------------
obj <- LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                spline = ~spl1D(x, nseg = 17, xlim = c(0,1)),
                data = dat, 
                family = multinomial())
summary(obj)

## ----makePredictions----------------------------------------------------------
sRows <- rowSums(multiNom)
fr <- multiNom/sRows
dat_fr <- data.frame(x, fr)

x0 <- seq(0, 1, by = 0.01)
newdat <- data.frame(x = x0)
pred <- predict(obj, newdata = newdat)
head(pred, 3)

## ----makePlot-----------------------------------------------------------------
library(tidyr)
colnames(pred) <- c("x", "category", "y")
prob_true <- t(sapply(X=x0, FUN = function(x) { nonlinear(x, mu)}))
colnames(prob_true) <- names(mu)
df_true <- data.frame(x = x0, prob_true)
prob_true_lf <- df_true %>% gather(key = "category",value="y", A:D)
dat_fr_lf <- dat_fr %>% gather(key = "category",value="y", A:D)
p1 <- ggplot(prob_true_lf, aes(x = x, y=y, color = category)) +
  geom_line(linetype='dashed') +
  geom_line(data=pred) +
  geom_point(data=dat_fr_lf)
p1


## ----oatsdata-----------------------------------------------------------------
## Load data.
data(oats.data)
head(oats.data, 3)

## ----define LVmodel-----------------------------------------------------------
## Add plot as factor.
oats.data$plotF <- as.factor(oats.data$plot)
## Define the precision matrix, see eqn (A2) in Boer et al (2020).
N <- nrow(oats.data)
cN <- c(1 / sqrt(N - 1), rep(0, N - 2), 1 / sqrt(N - 1))
D <- diff(diag(N), diff = 1)
Delta <- 0.5 * crossprod(D)
LVinv <- 0.5 * (2 * Delta + cN %*% t(cN))
## Add LVinv to list, with name corresponding to random term.
lGinv <- list(plotF = LVinv)

## ----modelLV------------------------------------------------------------------
obj7 <- LMMsolve(fixed = yield ~ rep + gen,
                 random = ~plotF, 
                 ginverse = lGinv, 
                 data = oats.data)

## ----summary_dev_VAR----------------------------------------------------------
round(deviance(obj7, relative = FALSE), 2)
summary(obj7, which = "variances")

## ----schmidtmodeloats---------------------------------------------------------
obj7b <- LMMsolve(fixed = yield ~ rep,
                 random = ~gen + rep:block, 
                 data = oats.data)

## ----effdim-------------------------------------------------------------------
EDdf <- effDim(obj7b)
EDdf

## ----genh2--------------------------------------------------------------------
round(subset(EDdf, Term == "gen")$Ratio, 3)

## ----APSIMdat-----------------------------------------------------------------
data(APSIMdat)
head(APSIMdat, 3)

## ----APSIMmodel---------------------------------------------------------------
obj8 <- LMMsolve(fixed = biomass ~ 1,
                 spline = ~spl1D(x = das, nseg = 50), 
                 data = APSIMdat)

## ----APSIMmodelSummary--------------------------------------------------------
summary(obj8)

## ----APSIMplot----------------------------------------------------------------
das_range <- range(APSIMdat$das)
newdat <- data.frame(das=seq(das_range[1], das_range[2], length = 300))
pred8 <- predict(obj8, newdata = newdat, se.fit = TRUE)
ggplot(data = APSIMdat, aes(x = das, y = biomass)) +
  geom_point(size = 1.0) +
  geom_line(data = pred8, aes(y = ypred), color = "blue", linewidth = 1) +
  geom_ribbon(data = pred8, aes(x = das,ymin = ypred-2*se, ymax = ypred+2*se),
              alpha = 0.2, inherit.aes = FALSE) + 
    labs(title = "APSIM biomass as function of time", 
       x = "days after sowing", y = "biomass (kg)") +
  theme_bw()

## ----APSIMDeriv---------------------------------------------------------------
plotDatDt <- predict(obj8, newdata = newdat, se.fit=TRUE, deriv="das")

ggplot(data = plotDatDt, aes(x = das, y = ypred)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(x=das, ymin = ypred-2*se, ymax = ypred+2*se),
              alpha=0.2, inherit.aes = FALSE) +  
  labs(title = "APSIM growth rate as function of time", 
       x = "days after sowing", y = "growth rate (kg/day)") +
  theme_bw()

## ----multipop-----------------------------------------------------------------
## Load data for multiparental population.
data(multipop)
head(multipop, 3)

## ----residualARG--------------------------------------------------------------
## Fit null model.
obj9 <- LMMsolve(fixed = pheno ~ cross, 
                 residual = ~cross, 
                 data = multipop)
dev0 <- deviance(obj9, relative = FALSE)

## ----groupOPTION--------------------------------------------------------------
## Fit alternative model - include QTL with probabilities defined in columns 3:5 
lGrp <- list(QTL = 3:5)
obj10 <- LMMsolve(fixed = pheno ~ cross, 
                 group = lGrp,
                 random = ~grp(QTL),
                 residual = ~cross,
                 data = multipop) 
dev1 <- deviance(obj10, relative = FALSE)

## ----approxPvalue-------------------------------------------------------------
## Deviance difference between null and alternative model.
dev <- dev0 - dev1
## Calculate approximate p-value. 
minlog10p <- -log10(0.5 * pchisq(dev, 1, lower.tail = FALSE))
round(minlog10p, 2)

## ----QTLeffects---------------------------------------------------------------
coef(obj10)$QTL

