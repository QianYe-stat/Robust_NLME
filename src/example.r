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

########################### data

dat <- read.table(here::here("data","hiv_dat.txt"), header=TRUE)

dat <- dat[complete.cases(dat),]
dat <- dat %>% 
  arrange(patid, day)%>%
  dplyr::select(patid, cohort, lgcopy, day) %>% 
  mutate(cohort=as.factor(cohort), day=day/max(day))



########################## Models 
dat1 <- groupedData(lgcopy~day|patid, data=dat)
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
start0 <- c(p1=exp(10),p2=exp(6),p3=5) # the starting value is based on the related lecture notes 

nls.fit  <- nls(lgcopy~nf(p1,p2,p3,day), data=dat1,start=start0)
start <- coef(nls.fit)
nlme.fit <- nlme(lgcopy~nf(p1,p2,p3,day),fixed = p1+p2+p3 ~1,random = p1+p3 ~1,
                 data =dat1,start=c(start))
summary(nlme.fit)
fixef(nlme.fit)
apply(ranef(nlme.fit),2,sd)


time <- unique(dat[,"day"])
time.fac <- as.factor(dat[,"day"])
lgsigma <- log(tapply(dat$lgcopy, time.fac, var))
sigma.fit <- lm(lgsigma~time)
coef(sigma.fit)[2]
########################## ROBUST NLME

# residual dispersion model:  
ndf <- 3
sigmaObject <- list(
  model=~1+day+(1+day|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=ndf,
  str.val=c(2*log(nlme.fit$sigma), coef(sigma.fit)[2]) 
)

# random effect dispersion model
ranCovObject <- list(
  varying.disp=~p3,
  ran.dist="inverse-Chi",
  df=ndf
)

# mean structure model:  
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nlmeObject <- list(
  nf = function(p1,p2,p3,t) p1+p2*exp(-p3*t),
  model= lgcopy ~ nf(p1,p2,p3,day),
  var=c("day"),
  fixed = p1+p2+p3 ~1,
  random = p1+p3 ~1,
  family='normal', 
  ran.dist='normal',
  sigma=sigmaObject,    # residual dispersion model (include residual random eff)
  ran.Cov=ranCovObject,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme.fit),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit),2,sd),  # starting value for fixed dispersion of random eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf, Inf) # upper bounds for  fixed dispersion of random eff
)





 
