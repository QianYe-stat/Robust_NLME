get_summary <- function(output.model, type){
  sd.out <- output.model$sd
  fix.out <- output.model$fixed
  
  na.index1 <- apply(sd.out,1,FUN=function(t){any(is.na(t))|any(t<0.001)})
  na.index2 <- apply(fix.out,1,FUN=function(t){any(t>=25)})
  
  na.index <- na.index1|na.index2
  
  cat("drop", sum(na.index), "invalid results")
   
  EST <- apply(output.model$fixed[!na.index,],2,mean)
  BIAS <- abs(EST-beta)
  SE.em <- apply(output.model$fixed[!na.index,],2,sd)
  SE <- apply(output.model$sd[!na.index,],2,FUN = function(t){sqrt(mean(t^2))})
  MSE <- apply(output.model$sqErr[!na.index,], 2, mean)
  Coverage <- apply(output.model$coverage[!na.index,],2,mean)
  
  
  EST.disp <- apply(output.model$dispersion[!na.index,], 2, mean)
  if(type=="nlme")  BIAS.disp <- abs(EST.disp-c(d,exp(alpha0)))
  if(type=="Rnlme") BIAS.disp <- abs(EST.disp-c(d,exp(alpha0), alpha1))
  
  SE.em.disp <- apply(output.model$dispersion[!na.index,], 2, sd)
  MSE.disp <- apply(output.model$dispersion.SqErr[!na.index,] , 2, mean)
  
  sum.model <- list(fixed=cbind(EST, BIAS,SE.em, SE, MSE, Coverage), 
                    dispersion=cbind(EST.disp, BIAS.disp, SE.em.disp, MSE.disp))
}

