##########################################
# The real data set contains cd4 data
# the model for sigma contains cd4* only
# cd4* is an unobserved true vaule, linked with a LME model
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
#nlme.fit1 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p2+p3 ~1,
#                 data =dat1,start=c(start))
nlme.fit1 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1+p3 ~1,
                  data =dat1,start=c(start))
anova(nlme.fit, nlme.fit1)
summary(nlme.fit1)

## model for cd4
cd4.fit <- lme(cd4~day+I(day^2), data=dat, random=~1|patid)
summary(cd4.fit)
fixef(cd4.fit)
as.numeric(VarCorr(cd4.fit)[,"StdDev"][1])

dat$cd4.pred <- fitted(cd4.fit)

fitted <- fitted(nlme.fit)
resid <- dat$lgcopy-fitted

time <- unique(dat[,"day"])
time.fac <- as.factor(dat[,"day"])
lgsigma <- log(tapply(resid, time.fac, var))
sigma.fit <- lm(lgsigma~time)
summary(sigma.fit)
anova(sigma.fit)
########################## Models for ROBUST NLME method: a_i ~ N(0,1)
# mean structure model:  
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2

sigma1 <- list(
  model=NULL,
  str.disp=as.numeric(VarCorr(cd4.fit)[,"StdDev"][2]),
  lower.disp=NULL,
  upper.disp=NULL,
  parName="xi"
)
lmeObject <- list(
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
  str.fixed=fixef(cd4.fit),  # starting value for fixed effect
  str.disp=as.numeric(VarCorr(cd4.fit)[,"StdDev"][1]),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf) # upper bounds for  fixed dispersion of random eff
  
)



# residual dispersion model:  
sigma2 <- list(
  model=~1+cd4.true+(1|patid),
  link='log',
  ran.dist="normal",
  str.fixed=c(2*log(nlme.fit$sigma), 0),
  str.disp=1,
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  dispName="siga",
  trueVal.model=list(var="cd4.true", model=lmeObject)
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
  sigma=sigma2,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme.fit),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObject1 <- list(
  nf = "nf1",
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1+p3 ~1,
  family='normal', 
  ran.dist='normal',
  fixName="beta",
  ranName="u",
  dispName="d",
  sigma=sigma2,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme.fit1),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit1),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects <- list(nlmeObject, lmeObject)
nlmeObjects1 <- list(nlmeObject1, lmeObject)

get_Jloglike(nlmeObjects1)

out <- Rnlme(nlmeObjects=nlmeObjects, long.data=dat, idVar="patid", sd.method="HL", 
             dispersion.SD = TRUE, independent.raneff = "byModel")
out1<- Rnlme(nlmeObjects=nlmeObjects1, long.data=dat, idVar="patid", sd.method="HL", 
             dispersion.SD = TRUE, independent.raneff = "byModel")

#out <- Rnlme(nlmeObjects=nlmeObjects, long.data=dat, idVar="patid", sd.method="HL", 
#             dispersion.SD = TRUE,  independent.raneff = TRUE)
out1$fixedest
out1$fixedSD
out1$dispersion
out1$dispSD
out1$AIC
out1$BIC

out1$AIC
out1$BIC

########################## Bootstrapping SE

group <- dat$patid  # grouping variable, e.g patient ID
uniqueID <- unique(group)   
n <- length(uniqueID)  # sample size
ni <- table(group)   # number of repeated measurements for each subject 
N <- nrow(dat)
t <- dat$day

# estimates from cd4.fit
gamma <- out$fixedest[c(6:8)]
xi <- out$dispersion["xi"]
sigma_b <- out$dispersion["sigb1"]

# estimates from Rnlme
d <- out$dispersion[c(1:3)]
Mat <- out$SIGMA
beta <- out$fixedest[c(1:3)]
alpha0 <- out$fixedest["alpha0"]
alpha1 <-  out$fixedest["alpha1"]
alpha <- c(alpha0, alpha1)
sigma_a <- out$dispersion["siga0"]
################ models
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2

sigma1.bt <- list(
  model=NULL,
  str.disp=xi,
  lower.disp=NULL,
  upper.disp=NULL,
  parName="xi"
)
lmeObject.bt <- list(
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
  sigma=sigma1.bt,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=gamma,  # starting value for fixed effect
  str.disp=c(sigma_b),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf) # upper bounds for  fixed dispersion of random eff
)

# residual dispersion model:  
sigma2.bt <- list(
  model=~1+cd4.true+(1|patid),
  link='log',
  ran.dist="normal",
  str.fixed=c(alpha0, alpha1),
  str.disp=sigma_a,
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  dispName="siga",
  trueVal.model=list(var="cd4.true", model=lmeObject.bt)
)

nlmeObject.bt <- list(
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
  sigma=sigma2.bt,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=beta,  # starting value for fixed effect
  str.disp=d,  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

nlmeObjects.bt <- list(nlmeObject.bt, lmeObject.bt)

## bootstrap
beta.est <- alpha.est <- gamma.est <- c()
disp.est <- c()
List.Rnlme <- NULL

rep <- 1
for(k in 1:rep){
  
  cat("This is run", k, "\n")
  
  Rnlme.fit  <-  0
  class(Rnlme.fit) <- "try-error"
  convg <- FALSE
  
  
  while(convg==FALSE | class(Rnlme.fit)=="try-error"){
    
    ## generate random effects
    a0 <- rnorm(n, sd=sigma_a)
    
    D <- diag(c(d, sigma_b)) %*% Mat %*% diag(c(d, sigma_b))
    ran <- rmvnorm(n, sigma=D)
    
    u <- ran[,c(1:3)]
    b1 <- ran[,4]
    
    
    simdat <- c()
    
    for(i in 1:n){
      indexi <- dat$patid==uniqueID[i] 
      nii <- ni[i]
      ti <- t[indexi]
      
      ## simulate CD4
      b1i <- b1[i]
      cd_errori <- rnorm(nii, sd=xi)
      cd_truei <- nf2(gamma[1]+b1i, gamma[2], gamma[3], ti)
      cd_obsi <- cd_truei+cd_errori
      
      ## get time-varying variance
      a0i <- a0[i]
      sdi <- sqrt(exp(alpha0+alpha1*cd_truei+a0i))
      errori <- rnorm(nii, sd=sdi)
      
      ## simulate lgcopy
      ui <- u[i,]
      tolEffi <- beta+ui
      lgcopyi <- nf1(tolEffi[1], tolEffi[2],tolEffi[3],ti)+errori
      
      
      dati <- data.frame(patid=uniqueID[i], day=ti, lgcopy=lgcopyi, cd4=cd_obsi)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    
    Rnlme.fit <- try(Rnlme(nlmeObject=nlmeObjects.bt, long.data=simdat, idVar="patid", 
                           independent.raneff = "byModel", iterMax=20))
    if(class(Rnlme.fit)!="try-error") convg <- Rnlme.fit$convergence
    }
  
  
  List.Rnlme[[k]] <- Rnlme.fit
  beta.est <- rbind(beta.est, Rnlme.fit$fixedest[c(1:3)])
  alpha.est <- rbind(alpha.est, Rnlme.fit$fixedest[c(4,5)])
  gamma.est <- rbind(gamma.est, Rnlme.fit$fixedest[c(6:8)])
  disp.est <- rbind(disp.est, Rnlme.fit$dispersion)
  
}

apply(beta.est, 2, sd)
apply(alpha.est, 2, sd)
apply(gamma.est, 2, sd)
apply(disp.est, 2, sd)

drop.index1.beta <- apply(beta.est,1,FUN=function(t){max(abs((t-beta)/beta))>0.15})
drop.index1.alpha <- apply(alpha.est,1,FUN=function(t){max(abs((t-alpha)/alpha))>0.15})
drop.index1 <- drop.index1.beta|drop.index1.alpha
rep-sum(drop.index1)

apply(beta.est[!drop.index1,], 2, sd)
apply(alpha.est[!drop.index1,], 2, sd)
apply(gamma.est[!drop.index1,], 2, sd)
apply(disp.est[!drop.index1,], 2, sd)

########################## latex table
library(xtable)
xtable(cbind(out$fixedest[c(1:3)], out$fixedSD[c(1:3)],
             apply(beta.est, 2, sd),
             apply(beta.est[!drop.index1,], 2, sd)
)
, type = "latex",digits = 3)

xtable(cbind(out$fixedest[c(4,5)], out$fixedSD[c(4,5)],
             apply(alpha.est, 2, sd),
             apply(alpha.est[!drop.index1,], 2, sd)
), type = "latex",digits = 3)

xtable(cbind(out$fixedest[c(6:8)], out$fixedSD[c(6:8)],
             apply(gamma.est, 2, sd),
             apply(gamma.est[!drop.index1,], 2, sd)
), type = "latex",digits = 3)

xtable(cbind(out$dispersion, out$dispSD,
             apply(disp.est, 2, sd),
             apply(disp.est[!drop.index1,], 2, sd)
), type = "latex",digits = 3)
