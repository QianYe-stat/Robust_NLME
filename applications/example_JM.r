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
summary(nlme.fit)

## model for cd4
cd4.fit <- lme(cd4~day+I(day^2), data=dat, random=~1|patid)
summary(cd4.fit)
fixef(cd4.fit)
as.numeric(VarCorr(cd4.fit)[,"StdDev"][1])

dat$cd4.pred <- fitted(cd4.fit)


########################## Models for ROBUST NLME method: a_i ~ N(0,1)
# mean structure model:  
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2

sigma1 <- list(
  model=NULL,
  str.fixed=as.numeric(VarCorr(cd4.fit)[,"StdDev"][2]),
  lower.fixed=NULL,
  upper.fixed=NULL,
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
  ran.dist="stdnormal",
  str.fixed=c(2*log(nlme.fit$sigma), 0),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
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


nlmeObjects <- list(nlmeObject, lmeObject)

out <- Rnlme(nlmeObjects=nlmeObjects, long.data=dat, idVar="patid", sd.method="HL", dispersion.SD = TRUE, iterMax=30)
out$fixedest
out$fixedSD
out$dispersion
out$dispSD

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
gamma <- 
xi <- 
sigma_b <- 

# estimates from Rnlme
d <- out$dispersion
Mat <- out$SIGMA
beta <- out$fixedest[c(1:3)]
alpha0 <- out$fixedest["alpha0"]
alpha1 <-  out$fixedest["alpha1"]
alpha <- c(alpha0, alpha1)
sigma_a <- 1
################ models
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2

sigma1.bt <- list(
  model=NULL,
  str.fixed=xi,
  lower.fixed=NULL,
  upper.fixed=NULL,
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
  ran.dist="stdnormal",
  str.fixed=c(alpha0, alpha1),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
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

