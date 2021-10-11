##########################################
# The real data set contains cd4 data
# the model for sigma contains cd4 only
#########################################
library(nlme)
library(tidyverse)
library(Deriv)
library(stringr)
library(LaplacesDemon)
library(purrr)
library(MASS)

rm(list=ls())
########################## source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)
set.seed(123)
########################### data

dat0 <- read.table(here::here("data","s315.txt"), header=TRUE)
names(dat0)

dat0 <- dat0[complete.cases(dat0),]
dat0 <- dat0 %>% 
  arrange(patid, day)%>%
  dplyr::select(patid, lgcopy, day, cd4) %>% 
  mutate(day=day/max(day), cd4=log(cd4))


ni <- tapply(dat0$day, dat0$patid, length)

sub_ID <- names(ni)[ni>=4]

dat <- subset(dat0, patid %in% sub_ID)

########################## Models for starting values

## mean structure model
dat1 <- groupedData(lgcopy~day|patid, data=dat)
plot(dat1)
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
start0 <- c(p1=10,p2=6,p3=5)
nls.fit  <- nls(lgcopy~nf(p1,p2,p3, day), data=dat1, start=start0)
start <- coef(nls.fit)

nlme.fit <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1+p2+p3 ~1,
                 data =dat1,start=c(start))
summary(nlme.fit)

## model for cd4
cd4.fit <- lme(cd4~day+I(day^2), data=dat, random=~1|patid)
summary(cd4.fit)

dat$cd4.pred <- fitted(cd4.fit)


########################## Models for ROBUST NLME method: a_i ~ N(0,1)

# residual dispersion model:  
sigmaObject2 <- list(
  model=~1+cd4.pred+(1|patid),
  link='log',
  ran.dist="stdnormal",
  str.fixed=c(2*log(nlme.fit$sigma), 0),
  lower.fixed=NULL,
  upper.fixed=NULL
)

# mean structure model:  
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)

nlmeObject2 <- list(
  nf = function(p1,p2,p3,t) p1+p2*exp(-p3*t),
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1+p2+p3 ~1,
  family='normal', 
  ran.dist='normal',
  sigma=sigmaObject2,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme.fit),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

get_nlme_loglike(nlmeObject2)


out.nor <- Rnlme(nlmeObject=nlmeObject2, long.data=dat, idVar="patid", sd.method="HL", dispersion.SD = TRUE)

out.nor$fixedest
out.nor$fixedSD
out.nor$SIGMA
out.nor$dispersion
out.nor$dispSD
out.nor$AIC
out.nor$BIC
out.nor$loglike_value

########################## Bootstrapping SE

group <- dat$patid  # grouping variable, e.g patient ID
uniqueID <- unique(group)   
n <- length(uniqueID)  # sample size
ni <- table(group)   # number of repeated measurements for each subject 
N <- nrow(dat)
t <- dat$day

# estimates from cd4.fit
cd4.coefs <- fixed.effects(cd4.fit)
cd4.sigma <- cd4.fit$sigma
cd4.randisp <- as.numeric(VarCorr(cd4.fit)[1,"StdDev"])

# estimates from Rnlme
d <- out.nor$dispersion[c(1:3)]
Mat <- out.nor$SIGMA
beta <- out.nor$fixedest
alpha0 <- out.nor$dispersion["alpha0"]
alpha1 <-  out.nor$dispersion["alpha1"]

############## models
sigmaObject.sim <- list(
  model=~1+cd4.pred+(1|patid),
  link='log',
  ran.dist="stdnormal",
  str.fixed=c(alpha0, alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL
)

nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)

nlmeObject.sim <- list(
  nf = function(p1,p2,p3,t) p1+p2*exp(-p3*t),
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1+p2+p3 ~1,
  family='normal', 
  ran.dist='normal',
  sigma=sigmaObject.sim,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)


beta.est.nor <- disp.est.nor <- gamma.est <- c()
List.Rnlme.nor <- NULL

rep <- 200
for(k in 1:rep){
  
  cat("This is run", k, "\n")
  
  Rnlme.fit <- lme.fit <-  0
  class(Rnlme.fit)<- class(lme.fit) <- "try-error"
  convg <- FALSE

  
  while(class(lme.fit)=="try-error"| convg==FALSE | class(Rnlme.fit)=="try-error"){
    
    cd4.ranef <- rnorm(n, sd=cd4.randisp)
    
    a0 <- rnorm(n, sd=1)
    D <- diag(d) %*% Mat %*% diag(d)
    u <- rmvnorm(n, sigma=D)
    
    simdat <- c()
    
    for(i in 1:n){
      indexi <- dat$patid==uniqueID[i] 
      nii <- ni[i]
      ti <- t[indexi]
      
      cd4.ranefi <- cd4.ranef[i]
      cd4.errori <- rnorm(nii, sd=cd4.sigma)
      cd4.truei <- cd4.coefs[1]+cd4.coefs[2]*ti+cd4.coefs[3]*ti^2+cd4.ranefi
      cd4.obsi <- cd4.truei+cd4.errori
      
      a0i <- a0[i]
      ui <- u[i,]
     
      sdi <- sqrt(exp(alpha0+alpha1*cd4.truei+a0i))
      errori <- rnorm(nii, sd=sdi)
      
      parami <- cbind(matrix(rep(beta+ui, nii), byrow=TRUE, ncol=length(beta)), ti)
      
      outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3], t[4])})
      
      dati <- data.frame(patid=uniqueID[i], lgcopy=outi+errori, day=ti, cd4=cd4.obsi)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    
    lme.fit <- try(lme(cd4~day+I(day^2), data=simdat, random=~1|patid))
    
    if(class(lme.fit)!="try-error"){
      simdat$cd4.pred <- fitted(lme.fit)
      
      Rnlme.fit <- try(Rnlme(nlmeObject=nlmeObject.sim, long.data=simdat, idVar="patid"))
      
      if(class(Rnlme.fit)!="try-error") convg <- Rnlme.fit$convergence
    }

  }
  
  List.Rnlme.nor[[k]] <- Rnlme.fit
  beta.est.nor <- rbind(beta.est.nor, Rnlme.fit$fixedest)
  disp.est.nor <- rbind(disp.est.nor, Rnlme.fit$dispersion[c(4,5)])
  gamma.est <- rbind(gamma.est, fixed.effects(lme.fit))
  
}
apply(beta.est.nor, 2, sd)
apply(disp.est.nor, 2, sd)
apply(gamma.est, 2, sd)

drop.index1.nor <- apply(beta.est.nor,1,FUN=function(t){max(abs((t-beta)/beta))>0.1})
rep-sum(drop.index1.nor)

apply(beta.est.nor[!drop.index1.nor,], 2, sd)
apply(disp.est.nor[!drop.index1.nor,], 2, sd)

drop.index2.nor <- apply(beta.est.nor,1,FUN=function(t){max(abs((t-beta)/beta))>0.05})
rep-sum(drop.index2.nor)

apply(beta.est.nor[!drop.index2.nor,], 2, sd)
apply(disp.est.nor[!drop.index2.nor,], 2, sd)

########################## latex table
library(xtable)
xtable(cbind(out.nor$fixedest, out.nor$fixedSD,
             apply(beta.est.nor, 2, sd),
             apply(beta.est.nor[!drop.index1.nor,], 2, sd),
             apply(beta.est.nor[!drop.index2.nor,], 2, sd)
), type = "latex",digits = 3)

xtable(cbind(out.nor$dispersion[c(4,5)], out.nor$dispSD[c(4,5)],
             apply(disp.est.nor, 2, sd),
             apply(disp.est.nor[!drop.index1.nor,], 2, sd),
             apply(disp.est.nor[!drop.index2.nor,], 2, sd)
), type = "latex",digits = 3)

xtable(cbind(cd4.coefs,
             summary(cd4.fit)$tTable[,"Std.Error"],
             apply(gamma.est, 2, sd),
             apply(gamma.est[!drop.index1.nor,], 2, sd),
             apply(gamma.est[!drop.index2.nor,], 2, sd)
), type = "latex",digits = 3)

save.image(here::here("results", "Rnlme_example3_bootstrap.RData"))
