rm(list=ls())
set.seed(123)
### read in the data
dat00 <- read_csv(here::here("data","raw_data.csv"))
names(dat00)
summary(dat00)
length(unique(dat00$ID))


### data process
dat0 <- dat00 %>% dplyr::select(ID, SEX, SEX_AT_BIRTH, FARVDT, DT, VLOAD, CD4, END_FOLLOW_DT, 
                       SUPP200_DT, SUPP200, REB200, REB200_DT, NON_BB_CLASSES_F, VD_BASE) %>% 
  mutate(SEX=as_factor(SEX),
         SEX_AT_BIRTH=as_factor(SEX_AT_BIRTH),
         NON_BB_CLASSES_F=as_factor(NON_BB_CLASSES_F),
         SUPP200=as.numeric(SUPP200),
         REB200=as.numeric(REB200),
         DT=as.Date(DT, format="%d-%b-%y"),
         FARVDT=as.Date(FARVDT, format="%d-%b-%y"),
         END_FOLLOW_DT=as.Date(END_FOLLOW_DT, format="%d-%b-%y"),
         SUPP200_DT=as.Date(SUPP200_DT, format="%d-%b-%y"),
         REB200_DT=as.Date(REB200_DT, format="%d-%b-%y"))%>% 
  rename(gender=SEX, sex=SEX_AT_BIRTH, ENDDT=END_FOLLOW_DT,
         supp=SUPP200, suppDT=SUPP200_DT, reb=REB200, rebDT=REB200_DT,
         patid=ID, cd4=CD4, trt=NON_BB_CLASSES_F) %>% 
  filter(DT>=FARVDT&DT<=ENDDT) %>% # only include records during the study period
  mutate(lgcopy=log10(VLOAD), cd4=sqrt(cd4), day=as.numeric(DT-FARVDT), time_to_SUPP=as.numeric(suppDT-FARVDT),
         FU_tol=ENDDT-FARVDT,
         supp=case_when(supp==99 ~ NA_real_,
                        TRUE ~ supp),
         reb=case_when(reb %in% c(88,99) ~ NA_real_,
                       reb %in% c(2,3) ~ 0,
                       TRUE ~ reb)) %>% 
  mutate(end_p1=case_when(
    supp==0 ~ ENDDT,
    is.na(supp) ~ ENDDT,
    supp==1 & reb==0 ~ ENDDT,
    supp==1 & reb==1 ~ rebDT,
    supp==1 & is.na(reb) ~ ENDDT
  ),
  end_p2=case_when(
    supp==0 ~ ENDDT,
    is.na(supp) ~ ENDDT,
    supp==1 ~ suppDT
  )) %>% 
  filter(day<=365) %>% # include records within the first year
  filter(DT <= end_p1) %>% # include records before rebound 
  filter(trt != "CAN/") %>% # exclude ppl with CAN regimen (only 5 ppl)
  # filter(FU_tol>=365.25) %>% # set min FU >=  1 year
  filter(VD_BASE>50) %>% # Naive patient with VL > 50
  mutate(trt=case_when( trt=="IIN/" ~ "IIN",
                        trt=="IIN/ENH/" ~ "IIN",
                        trt=="PI/" ~ "PI",
                        trt=="PI/ENH/" ~ "PI",
                        trt=="NNR/" ~ "NNR",
                        trt=="PI/IIN/ENH/" ~ "PI/IIN"
  )) %>% 
  mutate(day_raw=day, day=day/max(day)) %>% 
  mutate(trt_iin=if_else(trt=="IIN", 1, 0),
         trt_pi=if_else(trt=="PI", 1,0),
         trt_pi_iin=if_else(trt=="PI/IIN", 1,0))
  
uniqueID <- unique(dat0$patid)
length(uniqueID)

ni_vl <- tapply(dat0$lgcopy, dat0$patid, function(t){sum(!is.na(t))})
ni_cd  <-  tapply(dat0$cd4, dat0$patid, function(t){sum(!is.na(t))})
ni <-  pmin(ni_vl, ni_cd)
summary(ni)

long_ID <- names(ni)[ni>=4] # include patients with 4 or more measurements on lgcopy or cd4
dat <- dat0 %>% filter(patid %in% long_ID)

summary(dat)
table(dat$supp, exclude=NULL)
table(dat$reb, exclude=NULL)

uniqueID <- unique(dat$patid)
length(uniqueID)

ni_vl <- tapply(dat$lgcopy, dat$patid, function(t){sum(!is.na(t))})
ni_cd  <-  tapply(dat$cd4, dat$patid, function(t){sum(!is.na(t))})
ni <-  pmin(ni_vl, ni_cd)
summary(ni)
summary(ni_vl)
summary(ni_cd)
sum(!is.na(dat$lgcopy))
sum(!is.na(dat$cd4))
1-sum(!is.na(dat$cd4))/sum(!is.na(dat$lgcopy))
### Plot the trajectory
subID <- sample(uniqueID, 50)
subdat <- subset(dat, patid %in% subID) %>% filter(!is.na(lgcopy))

subdat1 <- groupedData(lgcopy~day_raw|patid, data=subdat)
  
plot(subdat1, xlab = "Time in days", ylab="Viral load in log10 scale")

ggplot(subdat, aes(day_raw, lgcopy, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('Time in days') + 
  ylab('Viral load in log10 scale')

ggplot(subdat, aes(day_raw, lgcopy, group=patid)) +geom_line()+
  geom_point(shape=1)+  theme(legend.position="bottom")+
  xlab('Time in days') + 
  ylab('Viral load in log10 scale')+
  facet_wrap(~trt, ncol=2)


### Plot the trajectory for NNR/
dat_nnr <- dat %>% filter(trt=="NNR/") %>% filter(!is.na(lgcopy))

add_dat <- tibble(times=seq(1:365)) %>% mutate(
  std_times=times/365,
  std_res=1.46+1.77*exp(-4.71*std_times),
  jm_res=1.54+2.22*exp(-10.03*std_times),
  patid=0)

colors <- c("Observed" = "black", "Standard" = "blue", "JM" = "red")
ggplot(dat_nnr, aes(day_raw, lgcopy, group=patid)) +
  geom_line()+
  geom_point(shape=1)+
  xlab('Time in days') + 
  ylab('Viral load in log10 scale') +
  geom_line(data=add_dat, mapping=aes(x=times, y=std_res, color="Standard"), size=1.2)+
  geom_line(data=add_dat, mapping=aes(x=times, y=jm_res, color="JM"), size=1.2)+
  labs(color="legend")+
  scale_color_manual(values = colors)

### Plot the trajectory for IIN/ENH
dat_iin_enh <- dat %>% filter(trt=="IIN/ENH/") %>% filter(!is.na(lgcopy))

add_dat <- tibble(times=seq(1:365)) %>% mutate(
  std_times=times/365,
  std_res=1.61+2.70*exp(-23.38*std_times),
  jm_res=1.54+2.67*exp(-31.33*std_times),
  patid=0)

colors <- c("Observed" = "black", "Standard" = "blue", "JM" = "red")
ggplot(dat_iin_enh, aes(day_raw, lgcopy, group=patid)) +
  geom_line()+
  geom_point(shape=1)+
  xlab('Time in days') + 
  ylab('Viral load in log10 scale') +
  geom_line(data=add_dat, mapping=aes(x=times, y=std_res, color="Standard"), size=1.2)+
  geom_line(data=add_dat, mapping=aes(x=times, y=jm_res, color="JM"), size=1.2)+
  labs(color="legend")+
  scale_color_manual(values = colors)



### save the output data
saveRDS(dat, here::here("data", "clean_data.rds"))
saveRDS(subID, here::here("data", "subID.rds"))
saveRDS(subdat, here::here("data", "subdat.rds"))
