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
library(matrixcalc)
library(lbfgs)
library(Matrix)
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
  dplyr::select(patid, lgcopy, day) %>% 
  mutate(day=day/max(day))


########################## Models for starting values
dat1 <- groupedData(lgcopy~day|patid, data=dat)
plot(dat1)
nf <- function(p1,p2,p3,t) p1+p2*exp(-p3*t)
start0 <- c(p1=exp(10),p2=exp(6),p3=5) # the starting value is based on the related lecture notes 

nls.fit  <- nls(lgcopy~nf(p1,p2,p3,day), data=dat1,start=start0)
start <- coef(nls.fit)
nlme.fit1 <- nlme(lgcopy~nf(p1,p2,p3,day),fixed = p1+p2+p3 ~1,random = p1+p2+p3 ~1,
                 data =dat1,start=c(start))
 
nlme.fit2 <- nlme(lgcopy~nf(p1,p2,p3,day),fixed = p1+p2+p3 ~1,random = p1+p3 ~1,
                  data =dat1,start=c(start))

time <- unique(dat[,"day"])
time.fac <- as.factor(dat[,"day"])
lgsigma <- log(tapply(dat$lgcopy, time.fac, var))
sigma.fit <- lm(lgsigma~time)

########################## Models for ROBUST NLME method

# residual dispersion model:  
ndf <- 5
sigmaObject1 <- list(
  model=~1+day+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=ndf,
  str.val=c(2*log(nlme.fit1$sigma), coef(sigma.fit)[2]),
  lower=NULL,
  upper=NULL
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
  str.fixed=fixef(nlme.fit1),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit1),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=NULL, # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

out1 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid")
out2 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid") # df=5
out3 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid")  # df=7
out4 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid")  # df=4

########################## Simulate one data set
group <- dat$patid  # grouping variable, e.g patient ID
uniqueID <- unique(group)   
n <- length(uniqueID)  # sample size
ni <- table(group)   # number of repeated measurements for each subject 
N <- nrow(dat)
t <- dat$day

d <- out2$dispersion[c(1:3)]
Mat <- out2$SIGMA
beta <- out2$fixedest
alpha0 <- out2$dispersion[4]
alpha1 <-  out2$dispersion[5]

## generate random effects
ndf <- 5
set.seed(123)
temp <- rchisq(n, df=ndf)
a0 <- log(ndf/temp)

temp <- rchisq(n, df=ndf)
expb3 <- ndf/temp

simdat <- c()
for(i in 1:n){
  a0i <- a0[i]
  expb3i <- expb3[i]
  nii <- ni[i]
  indexi <- dat$patid==uniqueID[i]
  ti <- t[indexi]

  sdi <- sqrt(exp(alpha0+alpha1*ti+a0i))
  errori <- rnorm(nii, sd=sdi)
  
  di <- d*c(1,1,sqrt(expb3i))
  
  Di <- diag(di) %*% Mat %*% diag(di)
  
  ui <- rmvnorm(1, sigma=Di)
  
  parami <- cbind(matrix(rep(beta+ui, nii), byrow=TRUE, ncol=3), ti)
  
  outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3], t[4])})
  
  dati <- tibble(patid=uniqueID[i], day=ti, lgcopy=outi+errori)
  
  simdat <- rbind(simdat, dati)
}

simdat <- simdat %>% 
  arrange(patid, day)

########################## plots: real VS simulated

gg2 <- ggplot(simdat, aes(day, lgcopy, group=patid)) +geom_line()+
  geom_point( shape=1)+  theme(legend.position="bottom")+
  xlab('time (re-scaled to [0,1])') + 
  ylab('lgcopy (Viral load in log10 scale)')+
  xlim(0,1)+
  ylim(1.8,6.8)+
  labs(title="Simulated data (k=5)")
print(gg2)

gg1 <- ggplot(dat, aes(day, lgcopy, group=patid)) +geom_line()+
  geom_point( shape=1)+  theme(legend.position="bottom")+
  xlab('time (re-scaled to [0,1])') + 
  ylab('lgcopy (Viral load in log10 scale)')+
  xlim(0,1)+
  ylim(1.8,6.8)+
  labs(title="Real data")
print(gg1)

ggarrange(gg1, gg2, ncol=2)

ggsave(here::here("figures", "ex1.pdf"))


#### 
 
plot_dat <- dat1 %>%  
  mutate(sim.lgcopy=simdat$lgcopy,
         patid=as.character(patid))

colors <- c("real" = "black", "simulated" = "blue")
ggplot(plot_dat, aes(day, lgcopy)) +
  geom_point(aes(day, sim.lgcopy, color="simulated"))+
  geom_line(aes(day, lgcopy, color="real"))+
  geom_line(aes(day, sim.lgcopy, color="simulated"))+
  geom_point()+ 
  xlab('time (re-scaled to [0,1])') + 
  ylab('lgcopy (Viral load in log10 scale)')+
  xlim(0,1)+
  ylim(1.8, 6.8)+
  facet_wrap(~patid, ncol=3)+
  labs(color="Legend")+
  scale_color_manual(values = colors)

ggsave(here::here("figures", "ex1_byID.pdf"))
