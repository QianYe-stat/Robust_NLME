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
dat0 <- read.table(here::here("data","s315.txt"), header=TRUE)
names(dat0)

dat0 <- dat0[complete.cases(dat0),]
dat0 <- dat0 %>% 
  arrange(patid, day)%>%
  dplyr::select(patid, lgcopy, day, cd4) %>% 
  mutate(day=day/max(day), cd4=log(cd4))

 
ni <- tapply(dat0$day, dat0$patid, length)

sub_ID <- names(ni)[ni>=6]

dat <- subset(dat0, patid %in% sub_ID)

########################## Models for starting values
dat1 <- groupedData(lgcopy~day|patid, data=dat)
plot(dat1)
#nf <- function(p1,p2,p3,p4, t, cd) p1+p2*exp(-(p4+p3*cd)*t)
nf <- function(p1,p2,p3, t, cd) p1+p2*t+p3*cd 
#start0 <- c(p1=exp(10),p2=exp(6),p3=1, p4=5) # the starting value is based on the related lecture notes 
start0 <- c(p1=exp(10),p2=exp(6),p3=1 ) # the starting value is based on the related lecture notes 
#nls.fit  <- nls(lgcopy~nf(p1,p2,p3, p4, day, cd4), data=dat1, start=start0)
nls.fit  <- nls(lgcopy~nf(p1,p2,p3, day, cd4), data=dat1, start=start0)
start <- coef(nls.fit)

fitted <- fitted(nls.fit)
resid <- dat$lgcopy-fitted

#nlme.fit <- nlme(lgcopy~nf(p1,p2,p3,p4, day, cd4),fixed = p1+p2+p3+p4 ~1,random = p1+p2+p4 ~1,
#                  data =dat1,start=c(start))

nlme.fit <- nlme(lgcopy~nf(p1,p2,p3, day, cd4),fixed = p1+p2+p3 ~1,random = p1+p2 ~1,
                  data =dat1,start=c(start))
summary(nlme.fit)


fitted <- fitted(nlme.fit)
resid <- dat$lgcopy-fitted

time <- unique(dat[,"day"])
time.fac <- as.factor(dat[,"day"])
lgsigma <- log(tapply(resid, time.fac, var))
sigma.fit <- lm(lgsigma~time)


# cd4
dat2 <- groupedData(cd4~day|patid, data=dat)
plot(dat2)
cd.fit <- lme(cd4~day+I(day^2), random=~1|patid, dat)
summary(cd.fit)
ranef(cd.fit)
########################## Models for ROBUST NLME method

# residual dispersion model:  
ndf <- 3
sigmaObject1 <- list(
  model=~1+day+(1|patid),
  link='log',
  ran.dist="inverse-Chi",
  df=ndf,
  str.val=c(2*log(nlme.fit$sigma), coef(sigma.fit)[2]) 
)

# random effect dispersion model
ranCovObject1 <- list(
  varying.disp=~p2,
  ran.dist="inverse-Chi",
  df=ndf
)


# mean structure model:  
#nf <- function(p1,p2,p3,p4, t, cd) p1+p2*exp(-(p4+p3*cd)*t)
nf <- function(p1,p2,p3, t, cd) p1+p2*t+p3*cd 
nlmeObject1 <- list(
  #nf = function(p1,p2,p3,p4, t, cd) p1+p2*exp(-(p4+p3*cd)*t),
  nf <- function(p1,p2,p3, t, cd) p1+p2*t+p3*cd,
  model= lgcopy ~ nf(p1,p2,p3, day, cd4),
  var=c("day", "cd4"),
  fixed = p1+p2+p3 ~1,
  random = p1+p2 ~1,
  family='normal', 
  ran.dist='normal',
  sigma=sigmaObject1,    # residual dispersion model (include residual random eff)
  ran.Cov=ranCovObject1,  # random effect dispersion model (include random random eff (double random eff))
  str.fixed=fixef(nlme.fit),  # starting value for fixed effect
  str.disp=apply(ranef(nlme.fit),2,sd),  # starting value for fixed dispersion of random eff
  lower.fixed=NULL, # lower bounds for fixed eff
  upper.fixed=rep(100,3), # upper bounds for fixed eff
  lower.disp=c(0,0,0), # lower bounds for fixed dispersion of random eff
  upper.disp=c(Inf,Inf,Inf) # upper bounds for  fixed dispersion of random eff
)

out1.ex2.df3 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid")
out1.ex2.df5 <- Rnlme(nlmeObject=nlmeObject1, long.data=dat, idVar="patid")

########################## Simulate one data set
group <- dat$patid  # grouping variable, e.g patient ID
uniqueID <- unique(group)   
n <- length(uniqueID)  # sample size
ni <- table(group)   # number of repeated measurements for each subject 
N <- nrow(dat)
t <- dat$day
cd <- dat$cd4

d <- out1.ex2.df3$dispersion[c(1:3)]
Mat <- out1.ex2.df3$SIGMA
beta <- out1.ex2.df3$fixedest
alpha0 <- out1.ex2.df3$dispersion[4]
alpha1 <-  out1.ex2.df3$dispersion[5]

## generate random effects
ndf <- 3
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
  cdi <- cd[indexi]
  
  sdi <- sqrt(exp(alpha0+alpha1*ti+a0i))
  errori <- rnorm(nii, sd=sdi)
  
  di <- d*c(1,1,sqrt(expb3i))
  
  Di <- diag(di) %*% Mat %*% diag(di)
  
  ui <- rmvnorm(1, sigma=Di)
  
  ui[4] <- ui[3]
  ui[3] <- 0
  
  parami <- cbind(matrix(rep(beta+ui, nii), byrow=TRUE, ncol=4), ti, cdi)
  
  outi <- apply(parami, 1, FUN=function(t){nf(t[1], t[2], t[3], t[4], t[5], t[6])})
  
  dati <- tibble(patid=uniqueID[i], day=ti, cd4=cdi, lgcopy=outi+errori)
  
  simdat <- rbind(simdat, dati)
}

simdat <- simdat %>% 
  arrange(patid, day)

########################## plots: real VS simulated
gg1 <- ggplot(dat, aes(day, lgcopy, group=patid)) +geom_line()+
  geom_point( shape=1)+  theme(legend.position="bottom")+
  xlab('time (re-scaled to [0,1])') + 
  ylab('lgcopy (Viral load in log10 scale)')+
  xlim(0,1)+
  ylim(1,8)+
  labs(title="Real data")
print(gg1)

gg2 <- ggplot(simdat, aes(day, lgcopy, group=patid)) +geom_line()+
  geom_point( shape=1)+  theme(legend.position="bottom")+
  xlab('time (re-scaled to [0,1])') + 
  ylab('lgcopy (Viral load in log10 scale)')+
  xlim(0,1)+
  ylim(1,8)+
  labs(title="Simulated data (k=3)")
print(gg2)

ggarrange(gg1, gg2, ncol=2)

ggsave(here::here("figures", "ex2.pdf"))


#### 

plot_dat <- dat %>%  
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
  ylim(1, 8)+
  facet_wrap(~patid, ncol=5)+
  labs(color="Legend")+
  scale_color_manual(values = colors)

ggsave(here::here("figures", "ex2_byID.pdf"), width=12, height=8)
