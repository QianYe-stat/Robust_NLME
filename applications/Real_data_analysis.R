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
########################### read in data

dat0 <- read.table(here::here("data","s315.txt"), header=TRUE)
names(dat0)


dat0 <- dat0 %>% 
  arrange(patid, day)%>%
  dplyr::select(patid, lgcopy, day, cd4) %>% 
  mutate(day1=day,day=day/max(day), cd4=log(cd4))
length(unique(dat0$patid))
table(dat0$day1)

#dat <- dat0[complete.cases(dat0),]
dat <- dat0

ni <- tapply(dat$day, dat$patid, length)

sub_ID <- names(ni)[ni>=3]

dat <- subset(dat, patid %in% sub_ID)
dat1 <- groupedData(lgcopy~day1|patid, data=dat)
########################### Descriptive 
summary(dat)
n=length(unique(dat$patid))
ni <- tapply(dat$day, dat$patid, length)

ID <- rep(c(1:n), ni)
dat$ID <- ID

summary(ni)

dat1 <- groupedData(lgcopy~day1|ID, data=dat)
plot(dat1, xlab = "Time in days", ylab="Viral load in log10 scale")

ggplot(dat, aes(day1, lgcopy, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('Time in days') + 
  ylab('Viral load in log10 scale')+
  xlim(0,92)

ggplot(dat, aes(day1, cd4, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('Time in days') + 
  ylab('CD4 in log scale')+
  xlim(0,92)
dat_cd <- groupedData(cd4~day1|ID, data=dat)

plot(dat_cd, xlab='Time in days', ylab='CD4 in log scale')
########################### Exploratory

### lme quadratic
fit_lm <- lmList(lgcopy~day+I(day^2)|patid, dat)
plot(intervals(fit_lm))

# Model 1: random intercept
lme1 <- lme(lgcopy~day+I(day^2), random=~1|patid, dat)
# Model 2: random intercept and linear term
lme2 <- lme(lgcopy~day+I(day^2), random=~1+day|patid, dat )
# Model 3: random intercept and quadratic term
lme3 <- lme(lgcopy~day+I(day^2), random=~1+I(day^2)|patid, dat )
# Model 4: random linear
lme4 <- lme(lgcopy~day+I(day^2), random=~0+day|patid, dat)
# Model 5: random quadratic term
lme5 <- lme(lgcopy~day+I(day^2), random=~0+I(day^2)|patid, dat )
# Model 6: random linear and quadratic term
lme6 <- lme(lgcopy~day+I(day^2), random=~0+day+I(day^2)|patid, dat )
# Model 7: random intercept, linear, and quadratic term
lme7 <- lme(lgcopy~day+I(day^2), random=~1+day+I(day^2)|patid, dat)

anova(lme1, lme2, lme3, lme4)
anova(lme1, lme2)
anova(lme1, lme3)
anova(lme1, lme4)
anova(lme4, lme5)
anova(lme2, lme4)
anova(lme3, lme4)

lme <- lme1
summary(lme)
plot(lme)
qqnorm(lme)
names(lme)
 
resid_lme <- residuals(lme, level = 1,type = "normalized")

ggplot(dat, aes(day, resid_lme, group=patid)) +geom_point()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('time (re-scaled to [0,1])') + 
  ylab('LME Residuals ')+
  xlim(0,1)

ggplot(dat, aes(day, resid_lme, group=patid)) +
  geom_point()+
  facet_wrap(~patid, ncol=9)

### nlme
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
start0 <- c(p1=10,p2=6,p3=5)
nls.fit  <- nls(lgcopy~nf(p1,p2,p3, day), data=dat, start=start0)
start <- coef(nls.fit)

# Model 1: random p1
nlme1 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1 ~1,
                 data =dat1,start=c(start))
# Model 2: random p1 p2
nlme2 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1+p2 ~1,
              data =dat1,start=c(start))
# Model 3: random p1 p3
nlme3 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1+p3 ~1,
              data =dat1,start=c(start))
# Model 4: random p1 p2 p3
nlme4 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1+p2+p3 ~1,
              data =dat1,start=c(start))
# Model 5: random p2
nlme5 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p2 ~1,
              data =dat1,start=c(start))
# Model 6: random p2 p3
nlme6 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p2+p3 ~1,
              data =dat1,start=c(start))
# Model 7: random p3
nlme7 <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p3 ~1,
              data =dat1,start=c(start))
anova(nlme1, nlme3)
anova(nlme1, nlme3)
anova(nlme3, nlme4)
anova(nlme3, nlme5)
anova(nlme3, nlme6)
anova(nlme3, nlme7)
summary(nlme3)

nlme <- nlme3
summary(nlme)
plot(nlme)

resid_nlme <- residuals(nlme, level = 1,type = "response")
ggplot(dat1, aes(day, resid_nlme, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('time (re-scaled to [0,1])') + 
  ylab('NLME Residuals')+
  xlim(0,1)

ggplot(dat1, aes(day, resid_nlme, group=patid)) +
  geom_point()+
  facet_wrap(~patid, ncol=9)

### Rnlme: use raw cd4 to model the residual variance
nf1 <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2
sigma <- list(
  model=~1+cd4+(1|patid),
  link='log',
  ran.dist="normal",
  str.fixed=c(2*log(nlme$sigma), 0),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  dispName="siga",
  str.disp=0.7
)


RnlmeObject <- list(
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
  sigma=sigma,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme),  # starting value for fixed effect
  str.disp=apply(ranef(nlme),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

RnlmeObjects=list(RnlmeObject)

dat_complete <- dat[complete.cases(dat),]
Rnlme.fit <- Rnlme(nlmeObjects=RnlmeObjects , long.data=dat_complete, 
      idVar="patid", sd.method="HL", dispersion.SD = TRUE,
      independent.raneff=FALSE)
Rnlme.fit$fixedest
Rnlme.fit$fixedSD
2*pnorm(abs(Rnlme.fit$fixedest/Rnlme.fit$fixedSD), lower.tail=FALSE)
Rnlme.fit$dispersion
Rnlme.fit$AIC

### two step Rnlme: use predicted cd4 to model the residual variance

# lme model for cd4
ggplot(dat_complete, aes(day1, cd4, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('Time in days') + 
  ylab('CD4 in log scale')+
  xlim(0,92)

cd4.fit1 <- lme(cd4~day, data=dat_complete, random=~1|patid, method = "ML")
cd4.fit2 <- lme(cd4~day+I(day^2), data=dat_complete, random=~1|patid, method="ML")
cd4.fit2 <- lme(cd4~day+I(day^2), data=dat_complete, random=~1|patid)
cd4.fit3 <- lme(cd4~day+I(day^2), data=dat_complete, random=~1+day|patid)
cd4.fit4 <- lme(cd4~day+I(day^2), data=dat_complete, random=~1+I(day^2)|patid)
cd4.fit5 <- lme(cd4~day+I(day^2), data=dat_complete, random=~1+day+I(day^2)|patid)

anova(cd4.fit1, cd4.fit2)
anova(cd4.fit2, cd4.fit3, cd4.fit3, cd4.fit5)

cd4.fit <- lme(cd4~day+I(day^2), data=dat_complete, random=~1|patid)

dat$cd4.pred <- predict(cd4.fit, newdata=dat)

# two step
sigma <- list(
  model=~1+cd4.pred+(1|patid),
  link='log',
  ran.dist="normal",
  str.fixed=c(2*log(nlme$sigma), 0),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  dispName="siga",
  str.disp=0.7
)


RnlmeObject <- list(
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
  sigma=sigma,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme),  # starting value for fixed effect
  str.disp=apply(ranef(nlme),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

Rnlme_TS=list(RnlmeObject)

Rnlme.TS.fit <- Rnlme(nlmeObjects=Rnlme_TS , long.data=dat, 
                   idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                   independent.raneff=FALSE, iterMax=20)
Rnlme.TS.fit$fixedest
Rnlme.TS.fit$fixedSD
2*pnorm(abs(Rnlme.TS.fit$fixedest/Rnlme.TS.fit$fixedSD), lower.tail=FALSE)
Rnlme.TS.fit$AIC
Rnlme.TS.fit$dispersion
summary(cd4.fit)

### JM with Rnlme: joint with a cd4 lme model, use true cd4 to model the residual variance
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
  str.fixed=c(2*log(nlme$sigma), 0),
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
  random = p1+p3 ~1,
  family='normal', 
  ran.dist='normal',
  fixName="beta",
  ranName="u",
  dispName="d",
  sigma=sigma2,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme),  # starting value for fixed effect
  str.disp=apply(ranef(nlme),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

Rnlme.JM <- list(nlmeObject, lmeObject)
get_Jloglike(Rnlme.JM)

Rnlme.JM.fit <- Rnlme(nlmeObjects=Rnlme.JM, long.data=dat_complete, idVar="patid", sd.method="HL", 
             dispersion.SD = TRUE, independent.raneff = "byModel")
Rnlme.JM.fit$fixedest
Rnlme.JM.fit$fixedSD
Rnlme.JM.fit$AIC

Rnlme.JM.fit.full <- Rnlme(nlmeObjects=Rnlme.JM, long.data=dat, idVar="patid", sd.method="HL", 
                      dispersion.SD = TRUE, independent.raneff = "byModel")
Rnlme.JM.fit.full$fixedest
Rnlme.JM.fit.full$fixedSD
2*pnorm(abs(Rnlme.JM.fit.full$fixedest/Rnlme.JM.fit.full$fixedSD), lower.tail=FALSE)
Rnlme.JM.fit.full$dispersion
Rnlme.JM.fit.full$AIC


get_sd_bootstrap1(Rnlme.fit=Rnlme.JM.fit.full, simdat=dat,at.rep=1 ,k.runs=2, big1=0.1, big2=0.15, 
                             independent.raneff = "byModel")

### JM with Rnlme: inverse-Chi
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
  ran.dist="inverse-Chi",
  str.fixed=c(2*log(nlme$sigma), 0),
  #str.disp=1,
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  df=3,
  #dispName="siga",
  trueVal.model=list(var="cd4.true", model=lmeObject)
)

nlmeObject <- list(
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
  str.fixed=fixef(nlme),  # starting value for fixed effect
  str.disp=apply(ranef(nlme),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

Rnlme.JM.invChi <- list(nlmeObject, lmeObject)


Rnlme.JM.fit.invChi <- Rnlme(nlmeObjects=Rnlme.JM.invChi, long.data=dat, idVar="patid", sd.method="HL", 
                      dispersion.SD = TRUE, independent.raneff = "byModel")

Rnlme.JM.fit.invChi$fixedest
Rnlme.JM.fit.invChi$fixedSD
save.image(here::here("results","Real_data_analysis.RData"))

round(unname(c(Rnlme.JM.fit.invChi$fixedest,
             Rnlme.JM.fit.invChi$dispersion[c("sigb1", "xi", "d1", "d3")])), 2)
saveRDS(Rnlme.JM.fit.full, here::here("results", "for real data analysis", "Rnlme.JM.fit.full.rds"))
saveRDS(Rnlme.JM.fit.invChi, here::here("results", "for real data analysis", "Rnlme.JM.fit.invChi.rds"))
saveRDS(dat, here::here("results", "for real data analysis", "realdata.rds"))

######################### Bootstrapping SE ###################
rm(list=ls())
########################## source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)
set.seed(123)
########################## read in models and data
dat <- readRDS(here::here("results", "for real data analysis", "realdata.rds"))
Rnlme.JM.fit.full <- readRDS(here::here("results", "for real data analysis", "Rnlme.JM.fit.full.rds"))
Rnlme.JM.fit.invChi <- readRDS(here::here("results", "for real data analysis", "Rnlme.JM.fit.invChi.rds"))

########################## get BT SE
# JM with Rnlme: joint with a cd4 lme model, use true cd4 to model the residual variance
JM.SE.bt <- get_sd_bootstrap1(Rnlme.fit=Rnlme.JM.fit.full, simdat=dat,at.rep=1 ,k.runs=2, big1=0.1, big2=0.15, 
                  independent.raneff = "byModel")
### JM with Rnlme: inverse-Chi
JM.invChi.SE.bt <- get_sd_bootstrap2(Rnlme.fit=Rnlme.JM.fit.invChi, simdat=dat,at.rep=1 ,k.runs=2, big1=0.1, big2=0.15, 
                                     independent.raneff = "byModel",df=3)

out <- list(JM.SE.bt=JM.SE.bt, JM.invChi.SE.bt=JM.invChi.SE.bt)
out
sebt <- readRDS(here::here("results","for real data analysis", "se.rds"))

(Rnlme.JM.fit.full$fixedSD.bt <- sebt$JM.SE.bt$se.bt)
2*pnorm(abs(Rnlme.JM.fit.full$fixedest/Rnlme.JM.fit.full$fixedSD.bt), lower.tail=FALSE)

(Rnlme.JM.fit.invChi$fixedSD.bt <- sebt$JM.invChi.SE.bt$se.bt)
2*pnorm(abs(Rnlme.JM.fit.invChi$fixedest/Rnlme.JM.fit.invChi$fixedSD.bt), lower.tail=FALSE)
