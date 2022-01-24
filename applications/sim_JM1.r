##########################################
# model setting can be found in /reports/Rnlme_example_3.pdf
# the model for sigma contains cd4 only
# for the Two-step (TS): the model for sigma contains cd4_pred
# for Joint models (JM): the model for sigma contains cd4* (unobserved true value)
# This simulation compare TS and JM
########################### load packages
library(nlme)
library(tidyverse)
library(Deriv)
library(stringr)
library(LaplacesDemon)
library(purrr)
library(MASS)
library(xtable)

rm(list=ls())
########################## source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)
######################### Simulation setting
rep <- 2
n <- 100
ni <- 15
N <- n*ni
ti <- seq(0, 1, length.out=ni)

patid <- rep(1:n, each=ni)
day <- rep(ti, n)
uniqueID <- seq(1:n)   

beta <- c(2.5, 3.0, 7.5)
#gamma <- c(5.2, 1.6, -1.2)
d <- c(0.4, 0.2, 1.0)
Mat <- matrix(c(1,0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), ncol=3)

alpha0 <- -8
alpha1 <-  1.5
sigma_a <- 1
#sigma_b <- 0.4
#xi <- 0.2

nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
#nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2

script <- "
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
"
writeLines(script, here::here("src","nf1.R"))


########################## model specification

####Two step

# residual dispersion model:  
sigmaObject <- list(
  model=~1+day+(1|patid),
  link='log',
  ran.dist="stdnormal",
  str.fixed=c(alpha0, alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a"
)


nlmeObject <- list(
  nf = "nf1",
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1+p2+p3 ~1,
  family='normal', 
  ran.dist='normal',
  fixName="beta",
  ranName="u",
  dispName="d",
  sigma=sigmaObject,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects=list(nlmeObject)


################################# Simulation runs
set.seed(123)
beta.est <- beta.sd <- c()
alpha.est <- alpha.sd <- c()

for(k in 1:rep){
  cat("This is run", k, "\n")
  
  Rnlme.fit <- 0
  class(Rnlme) <-  "try-error"
  convg <- FALSE
  
  while(class(Rnlme.fit)=="try-error" | convg==FALSE){
    
    ##########################  simulate data set
    
    ## generate random effects
    a0 <- rnorm(n, sd=sigma_a)
    
    
    D <- diag(d) %*% Mat %*% diag(d)
    u <- rmvnorm(n, sigma=D)
    
    simdat <- c()
    ## data set
    for(i in 1:n){
      
  
      ## get time-varying variance
      a0i <- a0[i]
      sdi <- sqrt(exp(alpha0+alpha1*ti+a0i))
      errori <- rnorm(ni, sd=sdi)
      
      ## simulate lgcopy
      ui <- u[i,]
      tolEffi <- beta+ui
      lgcopyi <- nf1(tolEffi[1], tolEffi[2],tolEffi[3],ti)+errori
      
      
      dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=lgcopyi)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    
    ########################## Modeling simdat
    Rnlme.fit <- try(Rnlme(nlmeObjects=nlmeObjects, long.data=simdat, 
                      idVar="patid", sd.method="HL", dispersion.SD = TRUE))

    if(class(Rnlme.fit)!="try-error") convg <- Rnlme.fit$convergence
      
  }
  ############### store output
  
  #### TS
  beta.est <- rbind(beta.est, Rnlme.fit$fixedest)
  beta.sd <- rbind(beta.sd, Rnlme.fit$fixedSD)
  
  alpha.index <- str_detect(names(Rnlme.fit$dispersion), "alpha")
  alpha.est <- rbind(alpha.est,Rnlme.fit$dispersion[alpha.index])
  alpha.sd <- rbind(alpha.sd, Rnlme.fit$dispSD[alpha.index])
  

  
  ##################### organize output 
  
  colnames(beta.est) <- colnames(beta.sd)  <- colnames(alpha.est) <- colnames(alpha.sd) <-  NULL
  
  out <- list(beta=list(True=beta,Est=beta.est, SD=beta.sd),
                 alpha=list(True=c(alpha0, alpha1),Est=alpha.est, SD=alpha.sd))

  
}
get_summary<- function(list){
  Est <- apply(list$Est, 2, mean)
  
  Bais_mat <- t(apply(list$Est, 1, FUN=function(t){t-list$True}))
  
  rBias <- abs(Est-list$True)/abs(list$True)*100
  
  rMSE <-  apply(Bais_mat, 2, FUN=function(t){mean(t^2)})/abs(list$True)*100
  
  SE.em <- apply(list$Est, 2, sd)
  
  SE <- apply(list$SD,2, FUN = function(t){sqrt(mean(t^2))})
  
  lower <- list$Est-1.96*list$SD
  upper <- list$Est+1.96*list$SD
  
  Coverage <- c()
  for(i in 1:length(list$True)){
    cov <- (lower[,i] <= list$True[i]) & (list$True[i]<= upper[,i])
    Coverage <- c(Coverage, mean(cov))
  }
  
  cbind(Est, rBias, rMSE, SE.em, SE, Coverage)
}

drop.index <- t(apply(out$beta$Est, 1, FUN=function(t){max(abs((t-out$beta$True)/out$beta$True))>0.05}))
sum(drop.index)
out1 <- out
out1$beta$Est <- out$beta$Est[!drop.index,]
out1$beta$SD <- out$beta$SD[!drop.index,]
out1$alpha$Est <- out$alpha$Est[!drop.index,]
out1$alpha$SD <- out$alpha$SD[!drop.index,]

res <- map(out1, get_summary)

res

save.image(here::here("results", "sim_JM1.RData"))


