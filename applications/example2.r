##########################################
# The real data set contains cd4 data
# the model for sigma contains both time and cd4
#########################################
library(nlme)
library(tidyverse)
library(Deriv)
library(stringr)
library(LaplacesDemon)
library(purrr)
library(MASS)

rm(list=ls())
#load(here::here("data", "Rnlme_bootstrap.RData"))
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
dat1 <- groupedData(lgcopy~day|patid, data=dat)
plot(dat1)
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
start0 <- c(p1=10,p2=6,p3=5)
nls.fit  <- nls(lgcopy~nf(p1,p2,p3, day), data=dat1, start=start0)
start <- coef(nls.fit)

nlme.fit <- nlme(lgcopy~nf(p1,p2,p3, day),fixed = p1+p2+p3 ~1,random = p1+p2+p3 ~1,
                  data =dat1,start=c(start))
summary(nlme.fit)
 

fitted <- fitted(nlme.fit)
resid <- dat$lgcopy-fitted

time <- unique(dat[,"day"])
time.fac <- as.factor(dat[,"day"])
lgsigma <- log(tapply(resid, time.fac, var))
sigma.fit <- lm(lgsigma~time)


########################## Models for ROBUST NLME method: a_i ~ inverse-Chisq

# residual dispersion model:  
ndf <- 5
sigmaObject1 <- list(
  model=~1+day+cd4+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=ndf,
  str.fixed=c(2*log(nlme.fit$sigma), coef(sigma.fit)[2], 0),
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
  str.fixed=fixef(nlme.fit),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

get_nlme_loglike(nlmeObject1)


out.df5 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid", sd.method="HL", dispersion.SD = TRUE)

# residual dispersion model:  
ndf <- 3
sigmaObject1 <- list(
  model=~1+day+cd4+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=ndf,
  str.fixed=c(2*log(nlme.fit$sigma), coef(sigma.fit)[2], 0),
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
  str.fixed=fixef(nlme.fit),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)



out.df3 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid", sd.method="HL", dispersion.SD = TRUE)

out.df5$fixedest
out.df5$fixedSD
out.df5$SIGMA
out.df5$dispersion
out.df5$dispSD
out.df5$AIC
out.df5$BIC
out.df5$loglike_value


out <- out.df3
out$fixedest
out$fixedSD
out$SIGMA
out$dispersion
out$dispSD
out$AIC
out$BIC
out$loglike_value

########################## Bootstrapping SE
group <- dat$patid  # grouping variable, e.g patient ID
uniqueID <- unique(group)   
n <- length(uniqueID)  # sample size
ni <- table(group)   # number of repeated measurements for each subject 
N <- nrow(dat)
t <- dat$day
cd <- dat$cd4

d <- out$dispersion[c(1:3)]
Mat <- out$SIGMA
beta <- out$fixedest
alpha0 <- out$dispersion["alpha0"]
alpha1 <-  out$dispersion["alpha1"]
alpha2 <- out$dispersion["alpha2"]

ndf <- 3
sigmaObject.sim <- list(
  model=~1+day+cd4+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=ndf,
  str.fixed=c(alpha0, alpha1, alpha2),
  lower.fixed=NULL,
  upper.fixed=NULL
)

# mean structure model:  
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


## generate random effects
ndf <- 3
beta.est <- disp.est <- c()
List.Rnlme <- NULL

rep <- 200
for(k in 1:rep){
  
  cat("This is run", k, "\n")
  
  Rnlme.fit <- 0
  class(Rnlme.fit)<- "try-error"
  convg <- FALSE
  rBias <- 1
  
  while(convg==FALSE | class(Rnlme.fit)=="try-error"){
  
    temp <- rchisq(n, df=ndf)
    a0 <- log(ndf/temp)
    D <- diag(d) %*% Mat %*% diag(d)
    u <- rmvnorm(n, sigma=D)
    
    simdat <- c()
    
    for(i in 1:n){
      a0i <- a0[i]
      ui <- u[i,]
      nii <- ni[i]
      indexi <- dat$patid==uniqueID[i]
      ti <- t[indexi]
      cdi <- cd[indexi]
      
      sdi <- sqrt(exp(alpha0+alpha1*ti+alpha2*cdi+a0i))
      errori <- rnorm(nii, sd=sdi)
      
      parami <- cbind(matrix(rep(beta+ui, nii), byrow=TRUE, ncol=length(beta)), ti)
      
      outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3], t[4])})
      
      dati <- data.frame(patid=uniqueID[i], lgcopy=outi+errori, day=ti, cd4=cdi)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    
    Rnlme.fit <- try(Rnlme(nlmeObject=nlmeObject.sim, long.data=simdat, idVar="patid"))
    if(class(Rnlme.fit)!="try-error") {
      convg <- Rnlme.fit$convergence
      rBias <- max(c(abs((Rnlme.fit$fixedest-beta)/beta)))
    }
  }
  
  List.Rnlme[[k]] <- Rnlme.fit
  #alpha0.ind <- names(Rnlme.fit$dispersion)=="alpha0"
  #Rnlme.fit$dispersion[alpha0.ind] <- exp(Rnlme.fit$dispersion[alpha0.ind])
  beta.est <- rbind(beta.est, Rnlme.fit$fixedest)
  disp.est <- rbind(disp.est, Rnlme.fit$dispersion)
  
}
apply(beta.est, 2, sd)
apply(disp.est, 2, sd)

drop.index1 <- apply(beta.est,1,FUN=function(t){max(abs((t-beta)/beta))>0.1})
rep-sum(drop.index1)

apply(beta.est[!drop.index1,], 2, sd)
apply(disp.est[!drop.index1,], 2, sd)

drop.index2 <- apply(beta.est,1,FUN=function(t){max(abs((t-beta)/beta))>0.05})
rep-sum(drop.index2)

apply(beta.est[!drop.index2,], 2, sd)
apply(disp.est[!drop.index2,], 2, sd)


########################## Models for ROBUST NLME method: a_i ~ N(0,1)

# residual dispersion model:  
sigmaObject2 <- list(
  model=~1+day+cd4+(1|patid),
  link='log',
  ran.dist="stdnormal",
  str.fixed=c(2*log(nlme.fit$sigma), coef(sigma.fit)[2], 0),
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

get_nlme_loglike(nlmeObject2)$randisp.df


out.nor <- Rnlme(nlmeObject=nlmeObject2, long.data=dat, idVar="patid", sd.method="HL", dispersion.SD = TRUE)

out.nor$fixedest
out.nor$fixedSD
out.nor$SIGMA
out.nor$dispersion
out.nor$AIC
out.nor$BIC
out.nor$loglike_value

############################ Booststraps ######################
group <- dat$patid  # grouping variable, e.g patient ID
uniqueID <- unique(group)   
n <- length(uniqueID)  # sample size
ni <- table(group)   # number of repeated measurements for each subject 
N <- nrow(dat)
t <- dat$day
cd <- dat$cd4

d <- out.nor$dispersion[c(1:3)]
Mat <- out.nor$SIGMA
beta <- out.nor$fixedest
alpha0 <- out.nor$dispersion["alpha0"]
alpha1 <-  out.nor$dispersion["alpha1"]
alpha2 <- out.nor$dispersion["alpha2"]

sigmaObject.sim <- list(
  model=~1+day+cd4+(1|patid),
  link='log',
  ran.dist="stdnormal",
  str.fixed=c(alpha0, alpha1, alpha2),
  lower.fixed=NULL,
  upper.fixed=NULL
)

# mean structure model:  
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


## generate random effects

beta.est.nor <- disp.est.nor <- c()
List.Rnlme.nor <- NULL

rep <- 200
for(k in 1:rep){
  
  cat("This is run", k, "\n")
  
  Rnlme.fit <- 0
  class(Rnlme.fit)<- "try-error"
  convg <- FALSE
  rBias <- 1
  
  while(convg==FALSE | class(Rnlme.fit)=="try-error"){
    
    a0 <- rnorm(n, sd=1)
    D <- diag(d) %*% Mat %*% diag(d)
    u <- rmvnorm(n, sigma=D)
    
    simdat <- c()
    
    for(i in 1:n){
      a0i <- a0[i]
      ui <- u[i,]
      nii <- ni[i]
      indexi <- dat$patid==uniqueID[i]
      ti <- t[indexi]
      cdi <- cd[indexi]
      
      sdi <- sqrt(exp(alpha0+alpha1*ti+alpha2*cdi+a0i))
      errori <- rnorm(nii, sd=sdi)
      
      parami <- cbind(matrix(rep(beta+ui, nii), byrow=TRUE, ncol=length(beta)), ti)
      
      outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3], t[4])})
      
      dati <- data.frame(patid=uniqueID[i], lgcopy=outi+errori, day=ti, cd4=cdi)
      
      simdat <- rbind(simdat, dati)
    }
    
    simdat <- simdat %>% arrange(patid, day)
    
    Rnlme.fit <- try(Rnlme(nlmeObject=nlmeObject.sim, long.data=simdat, idVar="patid"))
    if(class(Rnlme.fit)!="try-error") {
      convg <- Rnlme.fit$convergence
      rBias <- max(c(abs((Rnlme.fit$fixedest-beta)/beta)))
    }
  }
  
  List.Rnlme.nor[[k]] <- Rnlme.fit
  #alpha0.ind <- names(Rnlme.fit$dispersion)=="alpha0"
  #Rnlme.fit$dispersion[alpha0.ind] <- exp(Rnlme.fit$dispersion[alpha0.ind])
  beta.est.nor <- rbind(beta.est.nor, Rnlme.fit$fixedest)
  disp.est.nor <- rbind(disp.est.nor, Rnlme.fit$dispersion)
  
}
apply(beta.est.nor, 2, sd)
apply(disp.est.nor, 2, sd)

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
xtable(cbind(out.df5$fixedest, out.df5$fixedSD, 
             out.df3$fixedest, out.df3$fixedSD,
             apply(beta.est, 2, sd),
             apply(beta.est[!drop.index1,], 2, sd),
             apply(beta.est[!drop.index2,], 2, sd),
             out.nor$fixedest, out.nor$fixedSD,
             apply(beta.est.nor, 2, sd),
             apply(beta.est.nor[!drop.index1.nor,], 2, sd),
             apply(beta.est.nor[!drop.index2.nor,], 2, sd),
             fixed.effects(nlme.fit),
             summary(nlme.fit)$tTable[,"Std.Error"]
), type = "latex",digits = 3)
xtable(cbind(out.df5$dispersion, 
             out.df5$dispSD,
             out.df3$dispersion,
             out.df3$dispSD,
             apply(disp.est, 2, sd),
             apply(disp.est[!drop.index1,], 2, sd),
             apply(disp.est[!drop.index2,], 2, sd),
             out.nor$dispersion,
             out.nor$dispSD,
             apply(disp.est.nor, 2, sd),
             apply(disp.est.nor[!drop.index1.nor,], 2, sd),
             apply(disp.est.nor[!drop.index2.nor,], 2, sd),
             c(0.507, 0.332, 2.197,nlme.fit$sigma^2,1,1)
), type = "latex",digits = 3)

#save.image(here::here("data", "Rnlme_bootstrap.RData"))

########################## plots: real VS simulated
# plot_ID <- sample(names(ni),3)
# plot_dat <- subset(dat0, patid %in% plot_ID)
# 
# 
# gg1 <- ggplot(plot_dat, aes(day, lgcopy, group=patid)) +geom_line()+
#   geom_point( shape=1)+  theme(legend.position="bottom")+
#   xlab('time (re-scaled to [0,1])') + 
#   ylab('lgcopy (Viral load in log10 scale)')+
#   xlim(0,1)+
#   ylim(1,8)+
#   labs(title="Real data")
# print(gg1)
# ggsave(here::here("figures", "2.pdf"))
# 
# gg2 <- ggplot(simdat, aes(day, lgcopy, group=patid)) +geom_line()+
#   geom_point( shape=1)+  theme(legend.position="bottom")+
#   xlab('time (re-scaled to [0,1])') + 
#   ylab('lgcopy (Viral load in log10 scale)')+
#   xlim(0,1)+
#   ylim(1,8)+
#   labs(title="Simulated data (k=3)")
# print(gg2)
# 
# ggarrange(gg1, gg2, ncol=2)
# 
# ggsave(here::here("figures", "ex2.pdf"))
# 
# 
# #### 
# 
# plot_dat <- dat %>%  
#   mutate(sim.lgcopy=simdat$lgcopy,
#          patid=as.character(patid))
# 
# colors <- c("real" = "black", "simulated" = "blue")
# ggplot(plot_dat, aes(day, lgcopy)) +
#   geom_point(aes(day, sim.lgcopy, color="simulated"))+
#   geom_line(aes(day, lgcopy, color="real"))+
#   geom_line(aes(day, sim.lgcopy, color="simulated"))+
#   geom_point()+ 
#   xlab('time (re-scaled to [0,1])') + 
#   ylab('lgcopy (Viral load in log10 scale)')+
#   xlim(0,1)+
#   ylim(1, 8)+
#   facet_wrap(~patid, ncol=5)+
#   labs(color="Legend")+
#   scale_color_manual(values = colors)
# 
# ggsave(here::here("figures", "ex2_byID.pdf"), width=12, height=8)
