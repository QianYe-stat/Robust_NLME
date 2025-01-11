rm(list=ls())
dat <- readRDS(here::here("data","clean_data.rds"))
subID <- readRDS(here::here("data","subID.rds"))
subdat <-  readRDS(here::here("data","subdat.rds"))
########################## source all functions  
(file.sources = list.files(path=here::here("src"),pattern="*.R$"))
(file.sources <- paste0(here::here("src"), "/", file.sources))
sapply(file.sources,source)
set.seed(123)
####################################################
### nlme model
dat_vl <- dat %>% filter(!is.na(lgcopy)) 
dat1_vl <- groupedData(lgcopy~day|patid, data=dat_vl)

nf <- function(p1,p2,p3, b1, b2, b3, t, trt1, trt2, trt3 ) p1+p2*exp(-(p3+b1*trt1+b2*trt2+b3*trt3)*t)
start0 <- c(p1=10,p2=6,p3=5, b1=1, b2=2, b3=3)
nls.fit  <- nls(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ), data=dat_vl, start=start0)
start <- coef(nls.fit)

# Model 1: random p1
nlme1 <- nlme(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ),
              fixed = p1+p2+p3+b1+b2+b3 ~1,random = p1 ~1,
              data =dat1_vl,start=c(start))
# Model 2: random p1 p2
nlme2 <- nlme(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ),
              fixed = p1+p2+p3+b1+b2+b3 ~1,random = p1+p2 ~1,
              data =dat1_vl,start=c(start))
# Model 3: random p1 p3
nlme3 <- nlme(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ),
              fixed =  p1+p2+p3+b1+b2+b3 ~1,random = p1+p3 ~1,
              data =dat1_vl,start=c(start))
# Model 4: random p1 p2 p3
nlme4 <- nlme(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ),
              fixed =  p1+p2+p3+b1+b2+b3 ~1,random = p1+p2+p3 ~1,
              data =dat1_vl,start=c(start))
# Model 5: random p2
nlme5 <- nlme(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ),
              fixed =  p1+p2+p3+b1+b2+b3 ~1,random = p2 ~1,
              data =dat1_vl,start=c(start))
# Model 6: random p2 p3
nlme6 <- nlme(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ),
              fixed =  p1+p2+p3+b1+b2+b3 ~1,random = p2+p3 ~1,
              data =dat1_vl,start=c(start))
# Model 7: random p3
nlme7 <- nlme(lgcopy~nf(p1,p2,p3, b1, b2, b3, day, trt_iin, trt_pi, trt_pi_iin ),
              fixed =  p1+p2+p3+b1+b2+b3 ~1,random = p3 ~1,
              data =dat1_vl,start=c(start))

anova(nlme1, nlme2, nlme3, nlme6)
anova(nlme1, nlme2)
anova(nlme2, nlme3)
anova(nlme1, nlme3)
anova(nlme3, nlme6)
anova(nlme6, nlme7)

nlme <- nlme6
summary(nlme)
plot(nlme)

resid_nlme <- residuals(nlme, level = 1,type = "response")
sub_resid <- resid_nlme[names(resid_nlme) %in% subID]
ggplot(subdat, aes(day, sub_resid, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('time (re-scaled to [0,1])') + 
  ylab('NLME Residuals')+
  xlim(0,1)

ggplot(subdat, aes(day, sub_resid, group=patid)) +
  geom_point()+
  facet_wrap(~patid, ncol=9)


### two step Rnlme: use predicted cd4 to model the residual variance

# lme model for cd4
dat_cd <- dat %>% filter(!is.na(cd4)) 


subdat_cd <- subset(dat, patid %in% subID) %>% filter(!is.na(cd4))

ggplot(subdat_cd, aes(day_raw, cd4, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('Time in days') + 
  ylab('Square root of CD4')

cd4.fit1 <- lme(cd4~day, data=dat_cd, random=~1|patid)
cd4.fit2 <- lme(cd4~day+I(day^2), data=dat_cd, random=~1|patid, )
cd4.fit2 <- lme(cd4~day+I(day^2), data=dat_cd, random=~1|patid, )
cd4.fit3 <- lme(cd4~day+I(day^2), data=dat_cd, random=~1+day|patid)
cd4.fit4 <- lme(cd4~day+I(day^2), data=dat_cd, random=~1+I(day^2)|patid)
cd4.fit5 <- lme(cd4~day+I(day^2), data=dat_cd, random=~1+day+I(day^2)|patid)

anova(cd4.fit1, cd4.fit2)
anova(cd4.fit2, cd4.fit3, cd4.fit4, cd4.fit5)
anova(cd4.fit3, cd4.fit5)

cd4.fit <- cd4.fit5
summary(cd4.fit)
dat$cd4.pred <- predict(cd4.fit, newdata=dat)

# two step
nf1 <- function(p1,p2,p3, b1, b2, b3, t, trt1, trt2, trt3 ) p1+p2*exp(-(p3+b1*trt1+b2*trt2+b3*trt3)*t)
nf2 <- function(p1,p2,p3,t) p1+p2*t+p3*t^2
# sigma <- list(
#   model=~1+cd4.pred+(1|patid),
#   link='log',
#   ran.dist="normal",
#   str.fixed=c(2*log(nlme$sigma), 0),
#   lower.fixed=NULL,
#   upper.fixed=NULL,
#   fixName="alpha",
#   ranName="a",
#   dispName="siga",
#   str.disp=0.7
# )
sigma <- list(
  model=~1+cd4.pred+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  str.fixed=c(2*log(nlme$sigma), 0),
  lower.fixed=NULL,
  upper.fixed=NULL,
  fixName="alpha",
  ranName="a",
  #dispName="siga",
  #str.disp=0.7
  df=3
)


RnlmeObject <- list(
  nf = "nf1",
  model= lgcopy ~ nf(p1,p2,p3, b1, b2, b3, t, trt1, trt2, trt3 ),
  var=c("day", "trt_iin", "trt_pi", "trt_pi_iin"),
  fixed = p1+p2+p3+b1+b2+b3 ~1,
  random = p2+p3 ~1,
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
  upper.fixed=rep(100,6), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

Rnlme_TS=list(RnlmeObject)

Rnlme.TS.fit <- Rnlme(nlmeObjects=Rnlme_TS , long.data=as.data.frame(dat), 
                      idVar="patid", sd.method="HL", dispersion.SD = TRUE,
                      independent.raneff=FALSE, iterMax=50)
Rnlme.TS.fit$fixedest
Rnlme.TS.fit$fixedSD
2*pnorm(abs(Rnlme.TS.fit$fixedest/Rnlme.TS.fit$fixedSD), lower.tail=FALSE)
Rnlme.TS.fit$AIC
Rnlme.TS.fit$dispersion
summary(cd4.fit)

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
  random = p1+p2+p3 ~1,
  family='normal', 
  ran.dist='normal',
  fixName="gamma",
  ranName="b",
  dispName="sigb",
  sigma=sigma1,    # residual dispersion model (include residual random eff)
  ran.Cov=NULL,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(cd4.fit),  # starting value for fixed effect
  str.disp=as.numeric(VarCorr(cd4.fit)[,"StdDev"][c(1,2,3)]),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf, Inf, Inf) # upper bounds for  fixed dispersion of random eff
  
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
  model= lgcopy ~ nf(p1,p2,p3, b1, b2, b3, t, trt1, trt2, trt3 ),
  var=c("day", "trt_iin", "trt_pi", "trt_pi_iin"),
  fixed = p1+p2+p3+b1+b2+b3 ~1,
  random = p2+p3 ~1,
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
  upper.fixed=rep(100,6), # upper bounds for fixed eff
  lower.disp=c(0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

Rnlme.JM.invChi <- list(nlmeObject, lmeObject)


Rnlme.JM.fit.invChi <- Rnlme(nlmeObjects=Rnlme.JM.invChi, long.data=as.data.frame(dat), idVar="patid", sd.method="HL", 
                             dispersion.SD = TRUE, independent.raneff = "byModel", iterMax=50)

Rnlme.JM.fit.invChi$fixedest
Rnlme.JM.fit.invChi$fixedSD
2*pnorm(abs(Rnlme.JM.fit.invChi$fixedest/Rnlme.JM.fit.invChi$fixedSD), lower.tail=FALSE)
Rnlme.JM.fit.invChi$dispersion
save.image(here::here("results","02.1_model_fitting_adj_trt.RData"))
