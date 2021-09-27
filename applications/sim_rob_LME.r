<<<<<<< HEAD:applications/sim_rob_LME.r

rm(list=ls())
rep <- 200
n <- 100
ni <- 15
=======
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
rep <- 50
n <- 100
ni <- 10
>>>>>>> main:src/simulation2.r
N <- n*ni
ti <- seq(0, 1, length.out=ni)

patid <- rep(1:n, each=ni)
day <- rep(ti, n)

uniqueID <- seq(1:n)

<<<<<<< HEAD:applications/sim_rob_LME.r
beta <- c(6.0, -2.5)  
d <- c(0.42, 0.13)
Mat <- matrix(c(1,0.5, 0.5, 1), ncol=2)
=======
Mat <- matrix(c(1,0.1, 0.1, 0.1, 1, 0.1, 0.1, 0.1, 1), ncol=3)
>>>>>>> main:src/simulation2.r

alpha0 <- log(0.02)
alpha1 <-  3.4
ndf <- 3
nf <- function(p1,p2, t) p1+p2*t
script <- "
nf <- function(p1,p2, t) p1+p2*t
"
writeLines(script, here::here("src","nf.R"))
fit.df <- 3



########################## source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)


########################## model specification

# residual dispersion model:  
sigmaObject1 <- list(
  model=~1+day+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=fit.df,
  str.fixed=c(alpha0, alpha1) 
)


# mean structure model:  
nlmeObject1 <- list(
  nf = function(p1,p2, t) p1+p2*t ,
  model= lgcopy ~ nf(p1,p2, day),
  var=c("day"),
  fixed = p1+p2 ~1,
  random = p1+p2 ~1,
  family='normal', 
  ran.dist='normal',
  sigma=sigmaObject1,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,2), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

get_nlme_loglike(nlmeObject1)
set.seed(123)


beta.est.n <- beta.sd.n <- disp.est.n <- beta.COV.n <- beta.SqErr.n <- disp.SqErr.n <- c()
beta.est <- beta.sd <- disp.est <- beta.COV <- beta.SqErr<- disp.SqErr <- c()
<<<<<<< HEAD:applications/sim_rob_LME.r
=======
dat <- NULL 
## generate CD4
ai_cd <- rep(rnorm(n=n, sd=0.4), each=ni)
error_cd <- rnorm(N, sd=0.2)
cd4 <- 5.2+1.6*day-1.2*day^2+ai_cd+error_cd

>>>>>>> main:src/simulation2.r

for(k in 1:rep){
  cat("This is run", k, "\n")
  
  nlme.fit <- Rnlme.fit <- 0
  class(nlme.fit) <- class(Rnlme.fit)<- "try-error"
  convg <- FALSE
  
  while(class(nlme.fit)=="try-error" | convg==FALSE|class(Rnlme.fit)=="try-error"){
    ##########################  simulate data set
<<<<<<< HEAD:applications/sim_rob_LME.r
=======

>>>>>>> main:src/simulation2.r
    
    ## generate random effects
    temp <- rchisq(n, df=ndf)
    a0 <- log(ndf/temp) 
    
    D <- diag(d) %*% Mat %*% diag(d)
    u <- rmvnorm(n, sigma=D)

    
    simdat <- c()
    ## data set
    for(i in 1:n){
      
      a0i <- a0[i]
      ui <- u[i,]
      
      sdi <- sqrt(exp(alpha0+alpha1*ti+a0i))
      errori <- rnorm(ni, sd=sdi)
      
      parami <- cbind(matrix(rep((beta+ui), ni), byrow=TRUE, ncol=length(beta)), ti)
      
      outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3])})
      
      dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=outi+errori)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    #dat[[k]] <- simdat
    
    ########################## run nlme model
    simdat1 <- groupedData(lgcopy~day|patid, data=simdat)
<<<<<<< HEAD:applications/sim_rob_LME.r
    nf <-  function(p1,p2, t) p1+p2*t
    nlme.fit <- try(nlme(lgcopy~nf(p1,p2, day),fixed = p1+p2 ~1,random = p1+p2 ~1,
=======
    nf <- function(p1,p2,p3,p4, t, cd) p1+p2*exp(-(p4+p3*cd)*t)
    nlme.fit <- try(nlme(lgcopy~nf(p1,p2,p3,p4, day, cd4),fixed = p1+p2+p3+p4 ~1,random = p1+p2+p4 ~1,
>>>>>>> main:src/simulation2.r
                         data =simdat1,start=c(beta)))
    
    if(class(nlme.fit)!="try-error") {
      ########################## run Robust nlme model
      Rnlme.fit <- try(Rnlme(nlmeObject=nlmeObject1, long.data=simdat, idVar="patid", sd.method="aGH"))
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
  sigma.ind <- names(disp.nlme)=="Residual"
  disp.nlme[sigma.ind] <- (disp.nlme[sigma.ind])^2
  
  disp.SqErr.nlme <- (disp.nlme-c(d,exp(alpha0)))^2
  
  
  beta.est.n <- rbind(beta.est.n, beta.nlme)
  beta.sd.n <- rbind(beta.sd.n, sd.nlme) 
  beta.SqErr.n <- rbind(beta.SqErr.n, SqErr.nlme)
  beta.COV.n <- rbind(beta.COV.n, cover.nlme)
  
  
  disp.est.n <- rbind(disp.est.n, disp.nlme)
  disp.SqErr.n <- rbind(disp.SqErr.n, disp.SqErr.nlme)
  ################ output from Rnlme function
  alpha0.ind <- names(Rnlme.fit$dispersion)=="alpha0"
  Rnlme.fit$dispersion[alpha0.ind] <- exp(Rnlme.fit$dispersion[alpha0.ind])
  

  
  beta.lower <- Rnlme.fit$fixedest-1.96*Rnlme.fit$fixedSD
  beta.upper <- Rnlme.fit$fixedest+1.96*Rnlme.fit$fixedSD
  beta.cover <- (beta.lower<=beta) & (beta<=beta.upper)
<<<<<<< HEAD:applications/sim_rob_LME.r

=======
  cat("\n", beta.cover, "\n")
>>>>>>> main:src/simulation2.r
  
  beta.est <- rbind(beta.est, Rnlme.fit$fixedest)
  beta.sd <- rbind(beta.sd, Rnlme.fit$fixedSD)
  beta.SqErr <- rbind(beta.SqErr, (Rnlme.fit$fixedest-beta)^2)
  beta.COV <- rbind(beta.COV, beta.cover)
  
  disp.est <- rbind(disp.est, Rnlme.fit$dispersion)
  disp.SqErr <- rbind(disp.SqErr, (Rnlme.fit$dispersion-c(d,exp(alpha0),alpha1))^2)
  
}

output.nlme <- list(fixed=beta.est.n, sd=beta.sd.n, sqErr=beta.SqErr.n,
                    coverage=beta.COV.n, dispersion=disp.est.n, dispersion.SqErr=disp.SqErr.n)

output.Rnlme <- list(fixed=beta.est, sd=beta.sd, sqErr=beta.SqErr,
                     coverage=beta.COV, dispersion=disp.est, dispersion.SqErr=disp.SqErr)


get_summary <- function(output.model, type){
  sd.out <- output.model$sd
  fix.out <- output.model$fixed
  
  na.index1 <- apply(sd.out,1,FUN=function(t){any(is.na(t))|any(t<0.001)})
  na.index2 <- apply(fix.out,1,FUN=function(t){any(t>=25)})
  
  na.index <- na.index1|na.index2
  
  cat("drop", sum(na.index), "invalid results")
  
  EST <- apply(output.model$fixed[!na.index,],2,mean)
  BIAS <- abs(EST-beta)/abs(beta)*100
  SE.em <- apply(output.model$fixed[!na.index,],2,sd)
  SE <- apply(output.model$sd[!na.index,],2,FUN = function(t){sqrt(mean(t^2))})
  MSE <- apply(output.model$sqErr[!na.index,], 2, mean)/abs(beta)*100
  Coverage <- apply(output.model$coverage[!na.index,],2,mean)
  
  
  EST.disp <- apply(output.model$dispersion[!na.index,], 2, mean)
  if(type=="nlme") {
    BIAS.disp <- abs(EST.disp-c(d,exp(alpha0)))/abs(c(d,exp(alpha0)))*100
    MSE.disp <- apply(output.model$dispersion.SqErr[!na.index,] , 2, mean)/abs(c(d,exp(alpha0)))*100
  }
    
  if(type=="Rnlme") {
    BIAS.disp <- abs(EST.disp-c(d,exp(alpha0), alpha1))/abs(c(d,exp(alpha0), alpha1))*100
    MSE.disp <- apply(output.model$dispersion.SqErr[!na.index,] , 2, mean)/abs(c(d,exp(alpha0), alpha1))*100
  }
  
  SE.em.disp <- apply(output.model$dispersion[!na.index,], 2, sd)
  
  sum.model <- list(fixed=cbind(EST, BIAS,SE.em, SE, MSE, Coverage), 
                    dispersion=cbind(EST.disp, BIAS.disp, SE.em.disp, MSE.disp))
}


(sum.nlme <- get_summary(output.nlme, type='nlme'))
(sum.Rnlme <- get_summary(output.Rnlme, type="Rnlme"))


saveRDS(sum.nlme, here::here("data", "Rob_LME_aGH_nlmeRes.rds"))
saveRDS(sum.Rnlme, here::here("data", "Rob_LME_aGH_RnlmeRes.rds"))

xtable(cbind(sum.Rnlme$fixed, sum.nlme$fixed), type = "latex",digits = 3)
xtable(cbind(sum.Rnlme$dispersion, rbind(sum.nlme$dispersion,c(1,1,1,1))), type = "latex", digits = 3)


xtable(cbind(Rob_LME_aGH_RnlmeRes$fixed[,"EST"],
             Rob_LME_aGH_RnlmeRes$fixed[,"BIAS"],
             Rob_LME_aGH_RnlmeRes$fixed[,"MSE"], 
             Rob_LME_aGH_RnlmeRes$fixed[,"SE.em"],
             Rob_LME_aGH_RnlmeRes$fixed[,"SE"],
             Rob_LME_aGH_RnlmeRes$fixed[,"Coverage"],
             Rob_LME_aGH_nlmeRes$fixed[,"EST"],
             Rob_LME_aGH_nlmeRes$fixed[,"BIAS"],
             Rob_LME_aGH_nlmeRes$fixed[,"MSE"], 
             Rob_LME_aGH_nlmeRes$fixed[,"SE.em"],
             Rob_LME_aGH_nlmeRes$fixed[,"SE"],
             Rob_LME_aGH_nlmeRes$fixed[,"Coverage"]
), type = "latex",digits = 3)

<<<<<<< HEAD:applications/sim_rob_LME.r
xtable(cbind(Rob_LME_aGH_RnlmeRes$dispersion[,"EST.disp"],
             Rob_LME_aGH_RnlmeRes$dispersion[,"BIAS.disp"],
             Rob_LME_aGH_RnlmeRes$dispersion[,"MSE.disp"], 
             Rob_LME_aGH_RnlmeRes$dispersion[,"SE.em.disp"],
             Rob_LME_aGH_nlmeRes$dispersion[,"EST.disp"],
             Rob_LME_aGH_nlmeRes$dispersion[,"BIAS.disp"],
             Rob_LME_aGH_nlmeRes$dispersion[,"MSE.disp"], 
             Rob_LME_aGH_nlmeRes$dispersion[,"SE.em.disp"]
), type = "latex",digits = 3)
=======
>>>>>>> main:src/simulation2.r
