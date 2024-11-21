get_Jloglike <- function(nlmeObjects){
  
  k <- length(nlmeObjects) # number of lme/nlme models
  
  mu.loglike <- sigma.loglike <- ran.loglike <-  vector("list",k)
  
  Jfixed  <- Jdisp <- parSIGMA <- parNSIG <- str.fixed <- str.disp <- lower.fixed <- lower.disp <- upper.fixed <- upper.disp <- SIGMA.block <- c()
  
  for(i in 1:k){
    nlmeObject_i <- nlmeObjects[[i]]
    
    nlmeReturn <- get_nlme_loglike(nlmeObject=nlmeObject_i)
    
    lik <- nlmeReturn$loglike
    
    mu.loglike[i]<- lik$mu.loglike
    sigma.loglike[i] <- lik$sigma.loglike
    ran.loglike[i] <- lik$ran.loglike
    
    Jfixed <- c(Jfixed, nlmeReturn$fixed.par)
    #Jraneff <- c(Jraneff, nlmeReturn$ran.eff)
    Jdisp <- c(Jdisp, nlmeReturn$disp.par)
    
    parSIGMA <- c(parSIGMA, nlmeReturn$ran.eff[1:nlmeReturn$SIGMA.dim])
    parNSIG <- c(parNSIG, nlmeReturn$ran.eff[-(1:nlmeReturn$SIGMA.dim)])
    
    SIGMA.block <- c(SIGMA.block, nlmeReturn$SIGMA.dim)
    
    str.fixed <- c(str.fixed, nlmeReturn$str.fixed)
    str.disp <- c(str.disp,nlmeReturn$str.disp)
    
    lower.fixed <- c(lower.fixed, nlmeReturn$lower.fixed)
    lower.disp  <- c(lower.disp, nlmeReturn$lower.disp)
    
    upper.fixed <- c(upper.fixed, nlmeReturn$upper.fixed)
    upper.disp <- c(upper.disp,nlmeReturn$upper.disp)
  }
  
  Jraneff <- c(parSIGMA, parNSIG)
  
  #mu.loglike <- paste0("(", unlist(mu.loglike), ")", collapse="+")
  sigma.loglike <- paste0("(", unlist(sigma.loglike), ")", collapse="+")
  
  # likelihood of random effects
  SIGMA.dim <- length(parSIGMA)
  
  V.parSIGMA <- paste0("c(", paste0(parSIGMA, collapse = ","), ")")
  ran.loglike <- paste0("-0.5*",V.parSIGMA, "%*%invSIGMA%*%", V.parSIGMA,"+0.5*log(det(invSIGMA))-0.5*", SIGMA.dim,"*log(2*pi)")


  ## Joint likelihood
  Jloglike=list(mu.loglike=mu.loglike,sigma.loglike=sigma.loglike, ran.loglike=ran.loglike)
  

  result <- list(Jloglike=Jloglike, Jfixed=Jfixed, Jraneff=Jraneff, Jdisp=Jdisp,
       str.fixed=str.fixed, str.disp=str.disp, lower.fixed=lower.fixed, lower.disp=lower.disp,
       upper.fixed=upper.fixed, upper.disp=upper.disp,SIGMA.dim=SIGMA.dim, SIGMA.block=SIGMA.block)
  
  return(result)
}


 