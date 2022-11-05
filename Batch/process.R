## aggregate the results from sub-folders in VM
rm(list=ls())
library(xtable)
library(purrr)
res <- readRDS(here::here("Batch","sim15", "s1.rds"))
NLME.out <- res$NLME.out
TS.out <- res$TS.out
JM.out <- res$JM.out
alpha.NLME <- res$alpha.NLME

for(i in 2:10){
  filename <- paste0("s", i, ".rds")
  res <- readRDS(here::here("batch","sim15", filename))
  NLME.out <- map2(NLME.out, res$NLME.out, rbind)
  TS.out <- map2(TS.out, res$TS.out, rbind)
  JM.out <- map2(JM.out, res$JM.out, rbind)
  alpha.NLME <- c(alpha.NLME, res$alpha.NLME)
}
NLME.out$True <- NLME.out$True[1,]
TS.out$True <- TS.out$True[1,]
JM.out$True <- JM.out$True[1,]

big <- 15
rm.big <- function(out, big){
  out.bias <- t(apply(out$Est, 1, FUN=function(t){abs(t-out$True)/abs(out$True)*100}))
  out.rm <-  apply(out.bias,1,max)>big
  cat("\n",  deparse(substitute(out)) ,"output remove", sum(out.rm), "row with rBias >",big ,"%\n")
  out1 <- out
  out1$Est <- out$Est[!out.rm,]
  out1$SD <- out$SD[!out.rm,]
  if(!is.null(out$SD.BT)){
    out1$SD.BT <- out$SD.BT[!out.rm,]
    out1$SD.BT1 <- out$SD.BT1[!out.rm,]
    out1$SD.BT2 <- out$SD.BT2[!out.rm,]
  }
  
  return(list(out_df=out1, out.rm=out.rm))
}


get_summary<- function(out){
  Est <- apply(out$Est, 2, mean, na.rm=TRUE)
  
  Bais_mat <- t(apply(out$Est, 1, FUN=function(t){t-out$True}))
  
  rBias <- abs(Est-out$True)/abs(out$True)*100
  
  rMSE <-  apply(Bais_mat, 2, FUN=function(t){mean(t^2, na.rm=TRUE)})/abs(out$True)*100
  
  SE.em <- apply(out$Est, 2, sd, na.rm=TRUE)
  
  SE <- apply(out$SD,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
  
  
  lower <- out$Est-1.96*out$SD
  upper <- out$Est+1.96*out$SD
  
  Coverage <- c()
  for(i in 1:length(out$True)){
    cov <- (lower[,i] <= out$True[i]) & (out$True[i]<= upper[,i])
    Coverage <- c(Coverage, mean(cov, na.rm=TRUE))
  }
  
  if(!is.null(out$SD.BT)){
    SE.BT <- apply(out$SD.BT,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
    SE.BT1 <- apply(out$SD.BT1,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
    SE.BT2 <- apply(out$SD.BT2,2, FUN = function(t){sqrt(mean(t^2, na.rm=TRUE))})
    
    lower.bt <- out$Est-1.96*out$SD.BT
    upper.bt <- out$Est+1.96*out$SD.BT
    
    lower.bt1 <- out$Est-1.96*out$SD.BT1
    upper.bt1 <- out$Est+1.96*out$SD.BT1
    
    lower.bt2 <- out$Est-1.96*out$SD.BT2
    upper.bt2 <- out$Est+1.96*out$SD.BT2
    
    Coverage.bt <- Coverage.bt1 <-  Coverage.bt2 <- c()
    for(i in 1:length(out$True)){
      cov <- (lower.bt[,i] <= out$True[i]) & (out$True[i]<= upper.bt[,i])
      Coverage.bt <- c(Coverage.bt, mean(cov, na.rm=TRUE))
      
      cov <- (lower.bt1[,i] <= out$True[i]) & (out$True[i]<= upper.bt1[,i])
      Coverage.bt1 <- c(Coverage.bt1, mean(cov, na.rm=TRUE))
      
      cov <- (lower.bt2[,i] <= out$True[i]) & (out$True[i]<= upper.bt2[,i])
      Coverage.bt2 <- c(Coverage.bt2, mean(cov, na.rm=TRUE))
    }
    
  }
  
  res <- cbind(Est, rBias, rMSE, SE.em, SE, Coverage)
  
  if(!is.null(out$SD.BT)) res <- cbind(Est, rBias, rMSE, SE.em, SE, Coverage, 
                                       SE.BT, Coverage.bt,
                                       SE.BT1, Coverage.bt1,
                                       SE.BT2, Coverage.bt2)
  
  res
}


(nl <- get_summary(NLME.out))
mean(alpha.NLME)
(ts <- get_summary(TS.out))
(jm <- get_summary(JM.out))

NLME.out1 <- rm.big(NLME.out, big)$out_df
(nl1 <- get_summary(NLME.out1))
mean(alpha.NLME[!rm.big(NLME.out, big)$out.rm])


TS.out1 <- rm.big(TS.out, big)$out_df
(ts1 <- get_summary(TS.out1))

JM.out1 <- rm.big(JM.out, big)$out_df
(jm1 <- get_summary(JM.out1))

cat("\n xtable for original output \n ")
xtable(t(cbind(rbind(nl,0,0,0,0,0), ts, jm[,c(1:8)])), type = "latex",digits = 3)


cat("\n xtable for output with large rBias removed \n ")
xtable(cbind(rbind(nl1,0,0,0,0,0),ts1, jm1[,c(1:6,11,12)]), type = "latex",digits = 3)

runs.bt1 <- c(25.1,21.7,25.9,25,24,24.6, 23.4, 22.4, 21.9, 22.2)
runs.bt2 <- c(39.2,36.4,37.5,37.2,38.3,36.9,36.6, 36.1, 36.4, 35.7)
mean(runs.bt1)
mean(runs.bt2)
