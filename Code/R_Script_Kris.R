### This is an R script with functions that can be sourced into your R script.
### This can reduce the lines of code in your script
### Use the function source() to access this script. Run this line before using the pre-made functions
### For example: source("C:/Users/kreek001/OneDrive - Wageningen University & Research/R_Scripts")


## Show distribution, Q-Q plot Shapiro-Wilk test
Distr <- function(Independent.Variable) {
  par(mfrow = c(1 , 2))
  hist(Independent.Variable, prob = TRUE, col = "grey") #distribution histogram
  lines(density(Independent.Variable), col = "blue", lwd = 2)
  
  qqnorm(y = Independent.Variable, col = "blue", main = "Q-Q plot")
  qqline(y = Independent.Variable, col = "red")
  par(mfrow = c(1 , 1))
  print(shapiro.test(Independent.Variable)) #H0 = normal distribution assumed, when p < 0.05 H0 is rejected
}

### Test different distributions
library(fitdistrplus) #for testing distributions
library(actuar) #for testing distributions

Distr.No.Zero <- function(Ind.Var) {
  #norm, lnorm, exp, pois, cauchy, gamma, logis, nbinom, geom, beta, weibull from the stats package;
  n <- fitdist(Ind.Var, "norm")
  logn <- fitdist(log(Ind.Var), "norm")
  ln <- fitdist(Ind.Var, "lnorm") #no 0s in the data allowed
  ga <- fitdist(Ind.Var, "gamma") #no 0s in the data allowed
  logga <- fitdist(log(Ind.Var)+1.6, "gamma") #no 0s in the data allowed
  e <- fitdist(Ind.Var, "exp") 
  c <- fitdist(Ind.Var, "cauchy")
  li <- fitdist(Ind.Var, "logis")
  w <- fitdist(Ind.Var, "weibull") #no 0s allowed
  
  #invgamma, llogis, invweibull, pareto1, pareto from the actuar package.
  iga <- fitdist(Ind.Var, "invgamma") #observations need to be between 0 and 1
  lli <- fitdist(Ind.Var, "llogis") #no 0s allowed
  iw <- fitdist(Ind.Var, "invweibull") #no 0s allowed
  
  print(cbind(Fit.stats.Distr(n), Fit.stats.Distr(logn), Fit.stats.Distr(ln), Fit.stats.Distr(ga), Fit.stats.Distr(logga), 
              Fit.stats.Distr(e), Fit.stats.Distr(c), Fit.stats.Distr(li), Fit.stats.Distr(w), Fit.stats.Distr(iga), 
              Fit.stats.Distr(lli), Fit.stats.Distr(iw)))
  
  p.n <- plot(n, sub = "n")
  p.logn <- plot(logn, sub = "logn")
  p.ln <- plot(ln, sub = "ln")
  p.ga <- plot(ga, sub = "ga")
  p.logga <- plot(logga, sub = "logga") 
  p.e <- plot(e, sub = "e")
  p.c <- plot(c, sub = "c")
  p.li <- plot(li, sub = "li")
  p.w <- plot(w, sub = "w")
  p.iga <- plot(iga, sub = "iga") 
  p.lli <- plot(lli, sub = "lli")
  p.iw <- plot(iw, sub = "iw")
  
  distributions <- list(n, ln, ga, e, c, li, w, iga, lli, iw) # This does not work together with the log(normal) and log(gamma) data.
  plot.legend <- c("Normal", "lnorm", "gamma", "exponential", "cauchy", "llogis", 
                   "weibull", "inverse gamma", "llogis", "inverse weibull")
  
  cdfcomp(distributions, legendtext = plot.legend)
  denscomp(distributions, legendtext = plot.legend)
  qqcomp(distributions, legendtext = plot.legend)
  ppcomp(distributions, legendtext = plot.legend)
  gofstat(distributions)  
}


Distr.Zero <- function(Ind.Var) {
  ## Test several distributions
  #norm, lnorm, exp, pois, cauchy, gamma, logis, nbinom, geom, beta, weibull from the stats package;
  n <- fitdist(Ind.Var, "norm")
  logn <- fitdist(log(Ind.Var), "norm")
  ln <- fitdist(Ind.Var + 0.0001, "lnorm") #no 0s in the data allowed
  ga <- fitdist(Ind.Var + 0.0001, "gamma") #no 0s in the data allowed
  logga <- fitdist(log(Ind.Var)+1.6, "gamma") #no 0s in the data allowed
  e <- fitdist(Ind.Var, "exp") 
  c <- fitdist(Ind.Var, "cauchy")
  li <- fitdist(Ind.Var, "logis")
  w <- fitdist(Ind.Var + 0.0001, "weibull") #no 0s allowed
  
  #invgamma, llogis, invweibull, pareto1, pareto from the actuar package.
  iga <- fitdist(Ind.Var + 0.0001, "invgamma") #observations need to be between 0 and 1
  lli <- fitdist(Ind.Var + 0.0001, "llogis") #no 0s allowed
  iw <- fitdist(Ind.Var + 0.0001, "invweibull") #no 0s allowed
  
  print(cbind(Fit.stats.Distr(n), Fit.stats.Distr(logn), Fit.stats.Distr(ln), Fit.stats.Distr(ga), Fit.stats.Distr(logga), 
              Fit.stats.Distr(e), Fit.stats.Distr(c), Fit.stats.Distr(li), Fit.stats.Distr(w), Fit.stats.Distr(iga), 
              Fit.stats.Distr(lli), Fit.stats.Distr(iw)))
  
  p.n <- plot(n, sub = "n")
  p.logn <- plot(logn, sub = "logn")
  p.ln <- plot(ln, sub = "ln")
  p.ga <- plot(ga, sub = "ga")
  p.logga <- plot(logga, sub = "logga") 
  p.e <- plot(e, sub = "e")
  p.c <- plot(c, sub = "c")
  p.li <- plot(li, sub = "li")
  p.w <- plot(w, sub = "w")
  p.iga <- plot(iga, sub = "iga") 
  p.lli <- plot(lli, sub = "lli")
  p.iw <- plot(iw, sub = "iw")
  
  distributions <- list(n, ln, ga, e, c, li, w, iga, lli, iw) # This does not work together with the log(normal) and log(gamma) data.
  plot.legend <- c("Normal", "lnorm", "gamma", "exponential", "cauchy", "llogis", 
                   "weibull", "inverse gamma", "llogis", "inverse weibull")
  
  cdfcomp(distributions, legendtext = plot.legend)
  denscomp(distributions, legendtext = plot.legend)
  qqcomp(distributions, legendtext = plot.legend)
  ppcomp(distributions, legendtext = plot.legend)
  gofstat(distributions)  
  
}


## Comparing distributions ----
### Function to compare different distributions
Fit.stats.Distr <- function(Model){
  Model.Fit <- as.matrix(c(Model$aic, Model$bic, Model$loglik))
  rownames(Model.Fit) <- c("AIC", "BIC", "LogLik")
  colnames(Model.Fit) <- deparse(substitute(Model))
  return(Model.Fit)
}

## Comparing models ----
### Function by Daan to easily compare models for AIC, BIC and Log-likelihood of different models 
Fit.stats.Mod <- function(Model){
  Model.Fit <- as.matrix(c(AIC(Model), BIC(Model), logLik(Model)))
  rownames(Model.Fit) <- c("AIC", "BIC", "LogLik")
  colnames(Model.Fit) <- deparse(substitute(Model))
  return(Model.Fit)
}

### Compare model summaries of models with normal distribution
Fit.stats.Norm <- function(Model){
  #options(scipen = 999) # To make sure that there is not scientific notation
  summ <- summary(Model)
  Model.Fit <- as.matrix(c(AIC(Model), BIC(Model), logLik(Model), summ$fstatistic[[1]], summ$r.squared, summ$adj.r.squared))
  rownames(Model.Fit) <- c("AIC", "BIC", "LogLik", "F-stat", "R2", "Adj_R2")
  colnames(Model.Fit) <- deparse(substitute(Model))
  return(Model.Fit)
  #default_opts <- callr::r(function(){options("scipen")}); options(default_opts) #put scipen back to default (it does not work, why?)
}

### Compare model summaries of models with normal distribution
Fit.stats.Gam <- function(Model){
  #options(scipen = 999) # To make sure that there is not scientific notation
  summ <- summary(Model)
  Model.Fit <- as.matrix(c(AIC(Model), BIC(Model), logLik(Model), summ$deviance[[1]], summ$null.deviance, summ$adj.r.squared))
  rownames(Model.Fit) <- c("AIC", "BIC", "LogLik", "resdev", "nuldev")
  colnames(Model.Fit) <- deparse(substitute(Model))
  return(Model.Fit)
  #default_opts <- callr::r(function(){options("scipen")}); options(default_opts) #put scipen back to default (it does not work, why?)
}

#Model.Fit <- as.matrix(as.numeric(format(Model.Fit, scientific = FALSE)))


## Model assumptions ----
### Function to calculate residuals of a model
Res.sum1 <- function(Model) {
  fit.val  <- fitted(Model)
  res      <- residuals(Model) 
  p.res    <- residuals(Model, type = "pearson") #pearson residuals
  stun.res <- rstandard(Model) #stunderdised pearson residuals
  x <- cbind(fit.val, res, p.res, stun.res)
  print(summary(x))
  plot(p.res ~ fit.val, xlab="fitted values", ylab="Pearson residuals") # some values larger than 3 or smaller than -3, data may be over-dispersed
  abline(h = c(-3, -2, 0, 2, 3), col = c("red", "orange", "blue", "orange", "red"))
}

### Function to check model assumptions with the DHARMa package
library(DHARMa)
DHARMa.sum <- function(Model) {
  simulationOutput <- simulateResiduals(Model, n = 1000, seed = 256)
  par(mfrow = c(2, 2))
  cat("1. Kolmogorov-Smirnov test indicates deviations of the residuals from uniform distribution. H0: The model fits the data well; Ha: The model does not fit the data well")
  print(testUniformity(simulationOutput, plot = TRUE))
  plotResiduals(simulationOutput)               # Plot residuals
  cat("2. Test for over-dispersion. H0: No over-dispersion in the model; Ha: Over-dispersion in the model")
  print(testDispersion(simulationOutput, plot = TRUE)) # Tests over-dispersion: is observed data more/less dispersed than expected under fitted model
  cat("3. Test for outliers. H0: No ouliers in the model; Ha: Outliers in the model")
  print(testOutliers(simulationOutput, plot = TRUE))   # Tests for outliers
  par(mfrow = c(1, 1))
}
