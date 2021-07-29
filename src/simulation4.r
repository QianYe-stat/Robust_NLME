library(nlme)
library(tidyverse)
library(Deriv)
library(stringr)
library(LaplacesDemon)
library(purrr)
library(MASS)
library(mvtnorm)
library(tibble)
library(ggpubr)
library(Matrix)
rm(list=ls())
rep <- 100
n <- 100
ni <- 15
N <- n*ni
ti <- seq(0, 1, length.out=ni)

patid <- rep(1:n, each=ni)
day <- rep(ti, n)

uniqueID <- seq(1:n)   
d <- c(0.42, 0.13, 0.14)
sigma <- 0.2
Mat <- matrix(c(1,0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), ncol=3)
beta <- c(6.0, -2.5, -0.7)
alpha0 <- log(0.02)
alpha1 <-  3.4
ndf <- 3
nf <- function(p1,p2,p3, t, cd) p1+p2*t+p3*cd 
fit.df <- 3



########################## source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)


########################## model specification
sigmaObject1 <- list(
  model=NULL,
  str.val=sigma,
  lower=0,
  upper=Inf
)


# mean structure model:  
nlmeObject1 <- list(
  nf = function(p1,p2,p3, t, cd) p1+p2*t+p3*cd ,
  model= lgcopy ~ nf(p1,p2,p3, day, cd4),
  var=c("day", "cd4"),
  fixed = p1+p2+p3 ~1,
  random = p1+p2+p3 ~1,
  family='normal', 
  ran.dist='normal',
  sigma=sigmaObject1,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

get_nlme_loglike(nlmeObject1)
set.seed(123)

beta.est.n <- beta.sd.n <- disp.est.n <- beta.COV.n <- beta.SqErr.n <- disp.SqErr.n <- c()
beta.est <- beta.sd <- disp.est <- beta.COV <- beta.SqErr<- disp.SqErr <- c()
beta.sd1 <- beta.sd2 <- beta.sd3 <- beta.sd4 <- c()
beta.COV1 <- beta.COV2 <- beta.COV3 <- beta.COV4 <- c()
dat <- NULL 

## generate CD4
#ai_cd <- rep(rnorm(n=n, sd=0.4), each=ni)
#error_cd <- rnorm(N, sd=0.2)
#cd4 <- 5.2+1.6*day-1.2*day^2+ai_cd+error_cd

cd4 <- rnorm(N, mean=5.6, sd=0.5)

for(k in 1:rep){
  cat("This is run", k, "\n")
  
  nlme.fit <- Rnlme.fit <- 0
  class(nlme.fit) <- class(Rnlme.fit)<- "try-error"
  convg <- FALSE
  
  while(class(nlme.fit)=="try-error" | convg==FALSE|class(Rnlme.fit)=="try-error"){
    ##########################  simulate data set
    
    D <- diag(d) %*% Mat %*% diag(d)
    u <- rmvnorm(n, sigma=D)
    error <- rnorm(N, sd=sigma)
    
    simdat <- c()
    ## data set
    for(i in 1:n){
      ui <- u[i,]
      
      indexi <- patid==uniqueID[i]
      cdi <- cd4[indexi]
      errori <- error[indexi]
      
      parami <- cbind(matrix(rep((beta+ui), ni), byrow=TRUE, ncol=3), ti, cdi)
      
      outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3], t[4], t[5])})
      
      dati <- data.frame(patid=uniqueID[i], day=ti, cd4=cdi, lgcopy=outi+errori)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    #dat[[k]] <- simdat
    
    ########################## run nlme model
    simdat1 <- groupedData(lgcopy~day|patid, data=simdat)
    nf <-  function(p1,p2,p3, t, cd) p1+p2*t+p3*cd 
    nlme.fit <- try(nlme(lgcopy~nf(p1,p2,p3, day, cd4),fixed = p1+p2+p3 ~1,random = p1+p2+p3 ~1,
                         data =simdat1,start=c(beta)))
    
    if(class(nlme.fit)!="try-error") {
      ########################## run Robust nlme model
      Rnlme.fit <- try(Rnlme1(nlmeObject=nlmeObject1, long.data=simdat, idVar="patid"))
      if(class(Rnlme.fit)!="try-error") convg <- Rnlme.fit$convergence
    }
  }
  ################ output from nlme package
  nlme.summary <- summary(nlme.fit)
  names(nlme.fit)
  names(summary(nlme.fit))
  
  
  beta.nlme <- nlme.fit$coefficients$fixed
  sd.nlme <- nlme.summary$tTable[,"Std.Error"]
  SqErr.nlme <- (beta.nlme-beta)^2
  
  lower.nlme <- beta.nlme-1.96*sd.nlme
  upper.nlme <- beta.nlme+1.96*sd.nlme
  cover.nlme <- (lower.nlme <= beta) & (beta <= upper.nlme)
  
  disp.nlme <- as.numeric(VarCorr(nlme.fit)[,"StdDev"])
  names(disp.nlme) <- names(VarCorr(nlme.fit)[,"StdDev"])
 # sigma.ind <- names(disp.nlme)=="Residual"
 #  disp.nlme[sigma.ind] <- (disp.nlme[sigma.ind])^2
  
  disp.SqErr.nlme <- (disp.nlme-c(d,sigma))^2
  
  
  beta.est.n <- rbind(beta.est.n, beta.nlme)
  beta.sd.n <- rbind(beta.sd.n, sd.nlme) 
  beta.SqErr.n <- rbind(beta.SqErr.n, SqErr.nlme)
  beta.COV.n <- rbind(beta.COV.n, cover.nlme)
  
  
  disp.est.n <- rbind(disp.est.n, disp.nlme)
  disp.SqErr.n <- rbind(disp.SqErr.n, disp.SqErr.nlme)
  ################ output from Rnlme function
 # alpha0.ind <- names(Rnlme.fit$dispersion)=="alpha0"
#  Rnlme.fit$dispersion[alpha0.ind] <- exp(Rnlme.fit$dispersion[alpha0.ind])
  
  get_cov <- function(par,est,sd){
    lower <-  est-1.96*sd
    upper <-  est+1.96*sd
    cover <- (lower<=par) & (par<=upper)
  }
  
  #beta.lower <- Rnlme.fit$fixedest-1.96*Rnlme.fit$fixedSD
  #beta.upper <- Rnlme.fit$fixedest+1.96*Rnlme.fit$fixedSD
  #beta.cover <- (beta.lower<=beta) & (beta<=beta.upper)
  
  beta.cover1 <- get_cov(beta, Rnlme.fit$fixedest, Rnlme.fit$fixedSD1)
  beta.cover2 <- get_cov(beta, Rnlme.fit$fixedest, Rnlme.fit$fixedSD2)
  beta.cover3 <- get_cov(beta, Rnlme.fit$fixedest, Rnlme.fit$fixedSD3)
  beta.cover4 <- get_cov(beta, Rnlme.fit$fixedest, Rnlme.fit$fixedSD4)
  #cat("\n", beta.cover, "\n")
  
  beta.est <- rbind(beta.est, Rnlme.fit$fixedest)
  beta.sd1 <- rbind(beta.sd1, Rnlme.fit$fixedSD1)
  beta.sd2 <- rbind(beta.sd2, Rnlme.fit$fixedSD2)
  beta.sd3 <- rbind(beta.sd3, Rnlme.fit$fixedSD3)
  beta.sd4 <- rbind(beta.sd4, Rnlme.fit$fixedSD4)
  #beta.SqErr <- rbind(beta.SqErr, (Rnlme.fit$fixedest-beta)^2)
  beta.COV1 <- rbind(beta.COV1, beta.cover1)
  beta.COV2 <- rbind(beta.COV2, beta.cover2)
  beta.COV3 <- rbind(beta.COV3, beta.cover3)
  beta.COV4 <- rbind(beta.COV4, beta.cover4)
  
  disp.est <- rbind(disp.est, Rnlme.fit$dispersion)
  disp.SqErr <- rbind(disp.SqErr, (Rnlme.fit$dispersion-c(d,sigma))^2)
  
}
#output.nlme <- list(fixed=beta.est.n, sd=beta.sd.n, sqErr=beta.SqErr.n,
#                    coverage=beta.COV.n, dispersion=disp.est.n, dispersion.SqErr=disp.SqErr.n)

output.nlme <- list(fixed=beta.est.n, sd=beta.sd.n,
                    coverage=beta.COV.n, dispersion=disp.est.n )
#output.Rnlme <- list(fixed=beta.est, sd=beta.sd, sqErr=beta.SqErr,
 #                    coverage=beta.COV, dispersion=disp.est, dispersion.SqErr=disp.SqErr)

output.Rnlme <- list(fixed=beta.est, 
                     sd_TH=beta.sd1, coverage_TH=beta.COV1,
                     sd_TH_est=beta.sd2, coverage_TH_est=beta.COV2,
                     sd_HL=beta.sd3, coverage_HL=beta.COV3,
                     sd_HL_est=beta.sd4, coverage_HL_est=beta.COV4,
                     dispersion=disp.est)
map(output.nlme, function(t){apply(t,2,mean)})
map(output.Rnlme, function(t){apply(t,2,mean)})
saveRDS(map(output.nlme, function(t){apply(t,2,mean)}),  here::here("data", "nlmeCompare_n100_ni15_rep100.rds"))
saveRDS(map(output.Rnlme, function(t){apply(t,2,mean)}), here::here("data", "RnlmeCompare_n100_ni15_rep100.rds"))

(sum.nlme <- get_summary(output.nlme, type='nlme'))
(sum.Rnlme <- get_summary(output.Rnlme, type="Rnlme"))

saveRDS(sum.nlme, here::here("data", "nlmeRes_n=500"))
saveRDS(sum.Rnlme, here::here("data", "RnlmeRes_n100_linear.rds"))
