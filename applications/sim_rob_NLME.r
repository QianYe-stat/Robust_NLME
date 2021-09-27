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
rm(list=ls())
rep <- 500
n <- 100
ni <- 15
N <- n*ni
ti <- seq(0, 1, length.out=ni)

patid <- rep(1:n, each=ni)
day <- rep(ti, n)
uniqueID <- seq(1:n)   

beta <- c(3.4, 1.9, 15.7)
d <- c(0.6, 0.4, 2)
Mat <- matrix(c(1,0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), ncol=3)

alpha0 <- log(0.02)
alpha1 <-  5.2
ndf <- 5
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
script <- "
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
"
writeLines(script, here::here("src","nf.R"))
fit.df <- 5

########################## source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)


########################## model specification
# residual dispersion model:  
ndf <- 5
sigmaObject1 <- list(
  model=~1+day+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=fit.df,
  str.fixed=c(alpha0,alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL
)

# mean structure model:  
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)

nlmeObject1 <- list(
  nf = function(p1,p2,p3,t) p1+p2*exp(-p3*t),
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1+p2+p3 ~1,
  family='normal', 
  ran.dist='normal',
  sigma=sigmaObject1,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=NULL, # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)



set.seed(123)
beta.est.n <- beta.sd.n <- disp.est.n <- beta.COV.n <- beta.SqErr.n <- disp.SqErr.n <- c()
beta.est <- beta.sd.HL <- beta.sd.aGH <- disp.est <- beta.COV.HL <- beta.COV.aGH <- beta.SqErr<- disp.SqErr <- c()
List.Rnlme <- NULL
for(k in 1:rep){
  cat("This is run", k, "\n")
  
  nlme.fit <- Rnlme.fit <- 0
  class(nlme.fit) <- class(Rnlme.fit)<- "try-error"
  convg <- FALSE
#  rBias <- 1
  
  while(class(nlme.fit)=="try-error" | convg==FALSE|class(Rnlme.fit)=="try-error"){
    ##########################  simulate data set
    
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
      

      parami <- cbind(matrix(rep(beta+ui, ni), byrow=TRUE, ncol=length(beta)), ti)
      
      outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3], t[4])})
      
      dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=outi+errori)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
   # dat[[k]] <- simdat
    
    ########################## run nlme model
    simdat1 <- groupedData(lgcopy~day|patid, data=simdat)
    nf = function(p1,p2,p3,t) p1+p2*exp(-p3*t)
    nlme.fit <- try(nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1+p2+p3 ~1,
                         data =simdat1,start=c(beta)))
    
    if(class(nlme.fit)!="try-error") {
      ########################## run Robust nlme model
      Rnlme.fit <- try(Rnlme(nlmeObject=nlmeObject1, long.data=simdat, idVar="patid", sd.method="Both", ghsize=8))
      if(class(Rnlme.fit)!="try-error") {
        convg <- Rnlme.fit$convergence
     #   rBias <- max(c(abs((Rnlme.fit$fixedest-beta)/beta)))
      }
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
  List.Rnlme[[k]] <- Rnlme.fit
  alpha0.ind <- names(Rnlme.fit$dispersion)=="alpha0"
  Rnlme.fit$dispersion[alpha0.ind] <- exp(Rnlme.fit$dispersion[alpha0.ind])
  
  sd.HL <- Rnlme.fit$fixedSD$HL
  sd.aGH <- Rnlme.fit$fixedSD$aGH
  
  beta.lower.HL <- Rnlme.fit$fixedest-1.96*sd.HL
  beta.upper.HL <- Rnlme.fit$fixedest+1.96*sd.HL
  beta.cover.HL <- (beta.lower.HL<=beta) & (beta<=beta.upper.HL)
  
  beta.lower.aGH <- Rnlme.fit$fixedest-1.96*sd.aGH
  beta.upper.aGH <- Rnlme.fit$fixedest+1.96*sd.aGH
  beta.cover.aGH <- (beta.lower.aGH<=beta) & (beta<=beta.upper.aGH)
  
  beta.est <- rbind(beta.est, Rnlme.fit$fixedest)
  beta.SqErr <- rbind(beta.SqErr, (Rnlme.fit$fixedest-beta)^2)
  
  beta.sd.HL <- rbind(beta.sd.HL, sd.HL)
  beta.sd.aGH <- rbind(beta.sd.aGH, sd.aGH)
  
  beta.COV.HL <- rbind(beta.COV.HL, beta.cover.HL)
  beta.COV.aGH <- rbind(beta.COV.aGH, beta.cover.aGH)
  
  disp.est <- rbind(disp.est, Rnlme.fit$dispersion)
  disp.SqErr <- rbind(disp.SqErr, (Rnlme.fit$dispersion-c(d,exp(alpha0), alpha1))^2)
  
}
output.nlme <- list(fixed=beta.est.n, sd=beta.sd.n, sqErr=beta.SqErr.n,
                    coverage=beta.COV.n, dispersion=disp.est.n, dispersion.SqErr=disp.SqErr.n)
output.Rnlme <- list(fixed=beta.est, sd=list(HL=beta.sd.HL, aGH=beta.sd.aGH), sqErr=beta.SqErr,
                     coverage=list(HL=beta.COV.HL, aGH=beta.COV.aGH), 
                     dispersion=disp.est, dispersion.SqErr=disp.SqErr)


get_summary <- function(output.model, type, drop1=NULL, drop2=NULL, drop3=NULL){
  
  if(type=="nlme"){
    fix.out <- output.model$fixed
    drop.index <- apply(fix.out,1,FUN=function(t){max(abs((t-beta)/beta))>drop1})
  } 
  
  if(type=="Rnlme"){
    sd.out.aGH <- output.model$sd$aGH
    sd.out.HL <- output.model$sd$HL
    fix.out <- output.model$fixed
    
    drop.index1 <- apply(fix.out,1,FUN=function(t){max(abs((t-beta)/beta))>drop2})
    drop.index2 <- sd.out.aGH[,3]>drop3*sd.out.HL[,3]
    drop.index <- drop.index1 | drop.index2
  }
  
  
  cat("drop", sum(drop.index), "invalid results")
    
  EST <- apply(output.model$fixed[!drop.index,],2,mean)
  BIAS <- abs(EST-beta)/abs(beta)*100
  
  
  SE.em <- apply(output.model$fixed[!drop.index,],2,sd)
  MSE <- apply(output.model$sqErr[!drop.index,], 2, mean)/abs(beta)*100
  
  EST.disp <- apply(output.model$dispersion[!drop.index,], 2, mean)
  SE.em.disp <- apply(output.model$dispersion[!drop.index,], 2, sd)
  if(type=="nlme")  {
    SE <- apply(output.model$sd[!drop.index,],2,FUN = function(t){sqrt(mean(t^2))})
    Coverage <- apply(output.model$coverage[!drop.index,],2,mean)
    
    BIAS.disp <- abs(EST.disp-c(d,exp(alpha0)))/abs(c(d,exp(alpha0)))*100
    MSE.disp <- apply(output.model$dispersion.SqErr[!drop.index,], 2, mean)/abs(c(d,exp(alpha0)))*100
    
    sum.model <- list(fixed=cbind(EST, BIAS,MSE, SE.em,  SE, Coverage), 
                      dispersion=cbind(EST.disp, BIAS.disp, MSE.disp, SE.em.disp))
    
  }
  if(type=="Rnlme") {
    SE.HL <- apply(output.model$sd$HL[!drop.index,],2,FUN = function(t){sqrt(mean(t^2))})
    SE.aGH <- apply(output.model$sd$aGH[!drop.index,],2,FUN = function(t){sqrt(mean(t^2))})
    
    Coverage.HL <- apply(output.model$coverage$HL[!drop.index,],2,mean)
    Coverage.aGH <- apply(output.model$coverage$aGH[!drop.index,],2,mean)
    
    BIAS.disp <- abs(EST.disp-c(d,exp(alpha0), alpha1))/abs(c(d,exp(alpha0), alpha1))*100
    MSE.disp <- apply(output.model$dispersion.SqErr[!drop.index,], 2, mean)/abs(c(d,exp(alpha0), alpha1))*100 
    
    sum.model <- list(fixed=cbind(EST, BIAS,MSE, SE.em,  SE.HL , Coverage.HL, SE.aGH, Coverage.aGH), 
                    dispersion=cbind(EST.disp, BIAS.disp, MSE.disp, SE.em.disp))
  }
  sum.model
 
}


(sum.nlme <- get_summary(output.nlme, type='nlme', drop1=0.1))
(sum.Rnlme <- get_summary(output.Rnlme, type="Rnlme", drop2=0.05, drop3=5))

saveRDS(output.nlme, here::here("data", "Rob_NLME_nlme_output_8pts.rds"))
saveRDS(output.Rnlme, here::here("data", "Rob_NLME_Rnlme_output_8pts.rds"))
saveRDS(List.Rnlme, here::here("data", "Rob_NLME_Rnlme_list.rds"))

library(xtable)
xtable(cbind(sum.Rnlme$fixed, sum.nlme$fixed), type = "latex",digits = 3)
xtable(cbind(sum.Rnlme$dispersion, rbind(sum.nlme$dispersion,c(1,1,1,1))), type = "latex", digits = 3)


cbind(abs(output.Rnlme$fixed[,3]-15.7), output.Rnlme$sd$aGH[,3])
