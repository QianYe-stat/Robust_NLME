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
rep <- 5
n <- 100
ni <- 15
N <- n*ni
ti <- seq(0, 1, length.out=ni)

patid <- rep(1:n, each=ni)
day <- rep(ti, n)
uniqueID <- seq(1:n)   

beta <- c(2.5, 3.0, 7.5)
gamma <- c(5.2, 1.6, -1.2)
d <- 0.4
Mat <- matrix(c(1,0.5, 0.5, 1), ncol=2) # for d and b

alpha0 <- -8
alpha1 <-  1.5
sigma_a <- 1
sigma_b <- 0.4
xi <- 0.2

nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2

script <- "
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
"
writeLines(script, here::here("src","nf1.R"))

script <- "
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2
"
writeLines(script, here::here("src","nf2.R"))
########################## model specification

####Two step

# residual dispersion model:  
sigmaObject_TS <- list(
  model=~1+cd4.pred+(1|patid),
  link='log',
  ran.dist="stdnormal",
  str.fixed=c(alpha0, alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a"
)


nlmeObject_TS <- list(
  nf = "nf1",
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1 ~1,
  family='normal', 
  ran.dist='normal',
  fixName="beta",
  ranName="u",
  dispName="d",
  sigma=sigmaObject_TS,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects_TS=list(nlmeObject_TS)

#### JM

sigma1 <- list(
  model=NULL,
  str.fixed=xi,
  lower.fixed=NULL,
  upper.fixed=NULL,
  parName="xi"
)
lmeObject_JM <- list(
  nf = "nf2" ,
  model= cd4 ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1 ~1,
  family='normal', 
  ran.dist='normal',
  fixName="gamma",
  ranName="b",
  dispName="sigb",
  sigma=sigma1,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=gamma,  # starting value for fixed effect
  str.disp=c(sigma_b),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf) # upper bounds for  fixed dispersion of random eff
)

# residual dispersion model:  
sigma2 <- list(
  model=~1+cd4.true+(1|patid),
  link='log',
  ran.dist="stdnormal",
  str.fixed=c(alpha0, alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  trueVal.model=list(var="cd4.true", model=lmeObject_JM)
)

nlmeObject_JM <- list(
  nf = "nf1",
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1~1,
  family='normal', 
  ran.dist='normal',
  fixName="beta",
  ranName="u",
  dispName="d",
  sigma=sigma2,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects_JM <- list(nlmeObject_JM, lmeObject_JM)

################################# Simulation runs
set.seed(123)
beta.est.TS <- beta.sd.TS <-  beta.est.JM <- beta.sd.JM <- c()
alpha.est.TS <- alpha.sd.TS <- alpha.est.JM <- alpha.sd.JM <- c()
gamma.est.TS <- gamma.sd.TS <- gamma.est.JM <- gamma.sd.JM <- c()

for(k in 1:rep){
  cat("This is run", k, "\n")
  
  cd4.fit <- TS <- JM <- 0
  class(cd4.fit) <- class(TS) <- class(JM)<- "try-error"
  TSconvg <- JMconvg <-  FALSE
  
  while(class(cd4.fit)=="try-error" | class(JM)=="try-error" | class(TS)=="try-error"
        | TSconvg==FALSE| JMconvg==FALSE){
    
    ##########################  simulate data set
    
    ## generate random effects
    a0 <- rnorm(n, sd=sigma_a)
    
    D <- diag(c(d, sigma_b)) %*% Mat %*% diag(c(d, sigma_b))
    ran <- rmvnorm(n, sigma=D)
    
    u1 <- ran[,1]
    b1 <- ran[,2]
    
    simdat <- c()
    ## data set
    for(i in 1:n){
      
      ## simulate CD4
      b1i <- b1[i]
      cd_errori <- rnorm(ni, sd=xi)
      cdi_true <- nf2(gamma[1]+b1i, gamma[2], gamma[3], ti)
      cdi_obs <- cdi_true+cd_errori
      
      
      ## get time-varying variance
      a0i <- a0[i]
      sdi <- sqrt(exp(alpha0+alpha1*cdi_true+a0i))
      errori <- rnorm(ni, sd=sdi)
      
      ## simulate lgcopy
      ui <- c(u1[i],0,0)
      tolEffi <- beta+ui
      lgcopyi <- nf1(tolEffi[1], tolEffi[2],tolEffi[3],ti)+errori
      
      
      dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=lgcopyi, cd4=cdi_obs)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    
    ########################## Modeling simdat
    #### Two step
    
    cd4.fit <- try(lme(cd4~day+I(day^2), data=simdat, random=~1|patid))
    if(class(cd4.fit)!="try-error"){
      simdat$cd4.pred <- fitted(cd4.fit)
      
      cat("--Runing Two-step model\n")
      TS <- try(Rnlme(nlmeObjects=nlmeObjects_TS, long.data=simdat, 
                      idVar="patid", sd.method="HL", dispersion.SD = TRUE))
      cat("--done\n")
      if(class(TS)!="try-error") TSconvg <- TS$convergence
      
      if(class(TS)!="try-error" & TSconvg){
        #### JM
        cat("--Runing Joint model\n")
        JM <- try(Rnlme(nlmeObjects=nlmeObjects_JM, long.data=simdat, 
                        idVar="patid", sd.method="HL", dispersion.SD = TRUE, 
                        independent.raneff=TRUE))
        cat("--done\n")
        if(class(JM)!="try-error") JMconvg <- JM$convergence
      }
    }
  }
  ############### store output
  
  #### TS
  beta.index <- str_detect(names(TS$fixedest), "beta")
  beta.est.TS <- rbind(beta.est.TS, TS$fixedest[beta.index])
  beta.sd.TS <- rbind(beta.sd.TS, TS$fixedSD[beta.index])
  
  alpha.index <- str_detect(names(TS$fixedest), "alpha")
  alpha.est.TS <- rbind(alpha.est.TS,TS$fixedest[alpha.index])
  alpha.sd.TS <- rbind(alpha.sd.TS, TS$fixedSD[alpha.index])
  
  gamma.est.TS <- rbind(gamma.est.TS, fixef(cd4.fit))
  gamma.sd.TS <- rbind(gamma.sd.TS, summary(cd4.fit)$tTable[,"Std.Error"])
  
  #### JM
  
  beta.index <-  str_detect(names(JM$fixedest), "beta")
  beta.est.JM <- rbind(beta.est.JM, JM$fixedest[beta.index])
  beta.sd.JM <- rbind(beta.sd.JM, JM$fixedSD[beta.index])
  
  alpha.index <- str_detect(names(JM$fixedest), "alpha")
  alpha.est.JM <- rbind(alpha.est.JM,JM$fixedest[alpha.index])
  alpha.sd.JM <- rbind(alpha.sd.JM, JM$fixedSD[alpha.index])
  
  gamma.index <-  str_detect(names(JM$fixedest), "gamma")
  gamma.est.JM <- rbind(gamma.est.JM, JM$fixedest[gamma.index])
  gamma.sd.JM <-  rbind(gamma.sd.JM, JM$fixedSD[gamma.index])
  
  ##################### organize output 
  
  colnames(beta.est.TS) <- colnames(beta.sd.TS) <-  colnames(beta.est.JM) <- colnames(beta.sd.JM) <- colnames(alpha.est.TS) <- colnames(alpha.sd.TS) <- colnames(alpha.est.JM) <- colnames(alpha.sd.JM) <- colnames(gamma.est.TS) <- colnames(gamma.sd.TS) <- colnames(gamma.est.JM) <- colnames(gamma.sd.JM) <- NULL
  
  TS.out <- list(beta=list(True=beta,Est=beta.est.TS, SD=beta.sd.TS),
                 alpha=list(True=c(alpha0, alpha1),Est=alpha.est.TS, SD=alpha.sd.TS),
                 gamma=list(True=gamma,Est=gamma.est.TS, SD=gamma.sd.TS))
  JM.out <- list(beta=list(True=beta,Est=beta.est.JM, SD=beta.sd.JM),
                 alpha=list(True=c(alpha0, alpha1),Est=alpha.est.JM, SD=alpha.sd.JM),
                 gamma=list(True=gamma, Est=gamma.est.JM, SD=gamma.sd.JM))
  
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
map(TS.out, get_summary)
map(JM.out, get_summary)

#xtable(cbind(ts$beta, jm$beta), type = "latex",digits = 3)
#xtable(cbind(ts$alpha, jm$alpha), type = "latex",digits = 3)
#xtable(cbind(ts$gamma, jm$gamma), type = "latex",digits = 3)

save.image(here::here("results", "sim_JM.RData"))
#saveRDS(map(TS.out, get_summary), here::here("results", "TS_out.rds"))
#saveRDS(map(JM.out, get_summary), here::here("results", "JM_out.rds"))
