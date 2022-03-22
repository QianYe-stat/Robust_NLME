get_mu <- function(nlmeObject){
  ranCovObject <- nlmeObject$ranCovObject
  
  resp <- strsplit(as.character(nlmeObject$model), "~",  fixed=T)[[2]]
  resp <- str_trim(resp)
  
  rvX <- nlmeObject$var
  rvX <- str_trim(rvX)
  
  sp.fix <- strsplit(as.character(nlmeObject$fixed), "~",  fixed=T)[[2]]
  fix.comp <- strsplit(sp.fix, "+",  fixed=T)[[1]]
  fix.comp <- str_trim(fix.comp)
  
  sp.ran <- strsplit(as.character(nlmeObject$random), "~",  fixed=T)[[2]]
  ran.comp <- strsplit(sp.ran, "+",  fixed=T)[[1]]
  ran.comp <- str_trim(ran.comp)
  
  p <- length(fix.comp)  # dimension of fixed pars
  q <- length(ran.comp)  # dimension of random eff
  
  ran.ind <- c(1:p) * fix.comp %in% ran.comp
  
  if(!is.null(ranCovObject)){
    sp <- strsplit(as.character(ranCovObject$varying.disp), "~",  fixed=T)[[2]]
    randisp.comp <- strsplit(sp, "+",  fixed=T)[[1]]
    randisp.comp <- str_trim(randisp.comp) 
    randisp.ind <- c(1:p) * fix.comp %in% randisp.comp  # identify the random effects with varying dispersion
  }
  
  if(p>0){
    fixed <- paste0(nlmeObject$fixName, 1:p) # name fixed pars
  } else {fixed=NULL}
  
  if(q>0){
    raneff <- c()
    disp <- c()
    randisp <- c()
    
    for(i in 1:p){
      sub <- ran.ind[i]
      if(sub!=0){
        raneff <- c(raneff,paste0(nlmeObject$ranName, sub)) # name random effects
        disp <- c(disp, paste0(nlmeObject$dispName, sub))    # name the fixed dispersion of random eff
      } else {
        raneff <- c(raneff, 0)
        disp <- c(disp, 0)
      }
      
      rm(sub)
      if(!is.null(ranCovObject)){
        sub <- randisp.ind[i]
        if(sub!=0){
          randisp <- c(randisp, paste0(ranCovObject$ranName, sub)) # name the random effects for varying dispersion of random effect
        } else {                                  # double random effects
          randisp <- c(randisp, 0)
        }
      } else {
        randisp <- rep(0,q)
      }
    }
  } else {
    raneff=NULL
    disp=NULL
    randisp=NULL}
  
  disp.eff <- paste0(disp, rep("*", p), paste0("exp(0.5*", randisp, ")"))
  toteff <- paste0(fixed, rep("+", p), paste0(raneff, rep("*",p), disp.eff)) 
  
  mu <- paste0(nlmeObject$nf,"(", paste0(toteff, collapse = ","),"," ,paste0(rvX, collapse = ","), ")")
  
  return(mu)
  
}


get_nlme_loglike <- function(nlmeObject){
  
  ############## Return from get_info_sigma #####################
  if(is.null(nlmeObject$sigma$model)){
    sigma.raneff <- NULL  # random effect in residual dispersion (residual random effects)
    
    sigma.loglike <- NULL
    
    sigma.df <- NULL
    sigmaExpr <- paste0("(",nlmeObject$sigma$parName, "^2)" )
    
    sigma.fixed <- NULL 
    sigma.str <-  NULL
    sigma.lower <-  NULL
    sigma.upper <-  NULL
    
    sigma.disp <- nlmeObject$sigma$parName
    sigma.str.disp <- nlmeObject$sigma$str.disp
    sigma.lower.disp <- nlmeObject$sigma$lower.disp
    sigma.upper.disp <- nlmeObject$sigma$upper.disp
    
    if(is.null(sigma.lower.disp)) sigma.lower.disp <- 0
    if(is.null(sigma.upper.disp)) sigma.upper.disp <- Inf
    
  } else {
    sigmaInfo <- get_info_sigma(nlmeObject$sigma)
    sigma.raneff <- sigmaInfo$raneff  # random effect in residual dispersion (residual random effects)
    sigma.fixed <- sigmaInfo$fixed
    sigma.loglike <- sigmaInfo$loglike
    sigma.str <- sigmaInfo$str.val
    sigma.lower <- sigmaInfo$lower
    sigma.upper <- sigmaInfo$upper
    sigma.df <- sigmaInfo$df
    sigmaExpr <- sigmaInfo$sigmaExpr
    
    sigma.disp <- sigmaInfo$disp
    sigma.str.disp <- sigmaInfo$str.disp
    sigma.lower.disp <- sigmaInfo$lower.disp
    sigma.upper.disp <- sigmaInfo$upper.disp
  }
  
  ranCovObject <- nlmeObject$ran.Cov
  ##################
  
  resp <- strsplit(as.character(nlmeObject$model), "~",  fixed=T)[[2]]
  resp <- str_trim(resp)
  
  rvX <- nlmeObject$var
  rvX <- str_trim(rvX)
  
  sp.fix <- strsplit(as.character(nlmeObject$fixed), "~",  fixed=T)[[2]]
  fix.comp <- strsplit(sp.fix, "+",  fixed=T)[[1]]
  fix.comp <- str_trim(fix.comp)
  
  sp.ran <- strsplit(as.character(nlmeObject$random), "~",  fixed=T)[[2]]
  ran.comp <- strsplit(sp.ran, "+",  fixed=T)[[1]]
  ran.comp <- str_trim(ran.comp)
  
  p <- length(fix.comp)  # dimension of fixed pars
  q <- length(ran.comp)  # dimension of random eff
  
  ran.ind <- c(1:p) * fix.comp %in% ran.comp
  
  if(!is.null(ranCovObject)){
    sp <- strsplit(as.character(ranCovObject$varying.disp), "~",  fixed=T)[[2]]
    randisp.comp <- strsplit(sp, "+",  fixed=T)[[1]]
    randisp.comp <- str_trim(randisp.comp) 
    randisp.ind <- c(1:p) * fix.comp %in% randisp.comp  # identify the random effects with varying dispersion
  }
  
  if(p>0){
    fixed <- paste0(nlmeObject$fixName, 1:p) # name fixed pars
  } else {fixed=NULL}
  
  if(q>0){
    raneff <- c()
    disp <- c()
    randisp <- c()
    
    for(i in 1:p){
      sub <- ran.ind[i]
      if(sub!=0){
        raneff <- c(raneff,paste0(nlmeObject$ranName, sub)) # name random effects
        disp <- c(disp, paste0(nlmeObject$dispName, sub))    # name the fixed dispersion of random eff
      } else {
        raneff <- c(raneff, 0)
        disp <- c(disp, 0)
      }
      
      rm(sub)
      if(!is.null(ranCovObject)){
        sub <- randisp.ind[i]
        if(sub!=0){
          randisp <- c(randisp, paste0(ranCovObject$ranName, sub)) # name the random effects for varying dispersion of random effect
        } else {                                  # double random effects
          randisp <- c(randisp, 0)
        }
      } else {
        randisp <- rep(0,q)
      }
    }
  } else {
    raneff=NULL
    disp=NULL
    randisp=NULL}
  
  disp.eff <- paste0(disp, rep("*", p), paste0("exp(0.5*", randisp, ")"))
  toteff <- paste0(fixed, rep("+", p), paste0(raneff, rep("*",p), disp.eff)) 
  
  mu <- paste0(nlmeObject$nf,"(", paste0(toteff, collapse = ","),"," ,paste0(rvX, collapse = ","), ")")
  
  if (nlmeObject$family=="normal"){ 
    
    sigma <-  sigmaExpr
    mu.loglike <- paste0("-0.5*(", resp, "-", mu, 
                         ")^2/",sigma, "-0.5*log(",sigma,")-0.5*log(2*pi)")
    
  }
  
  raneff <- str_subset(raneff, "[^0]")
  disp <- str_subset(disp, "[^0]")
  randisp <- str_subset(randisp, "[^0]")
  
  if(!is.null(ranCovObject)){
    if(ranCovObject$ran.dist=="inverse-Chi"){
      randisp.loglike <- make_loglike_invChi(randisp, ranCovObject$df)
    }
  } 
  
  #####
    str.fixed <- c(nlmeObject$str.fixed, sigma.str)
    names(str.fixed) <- c(fixed, sigma.fixed)
    str.disp <- c(nlmeObject$str.disp, sigma.str.disp)
    names(str.disp) <- c(disp,sigma.disp)
  
  # if(q==1){
  #   ran.loglike <- make_loglike_normal(raneff, mean=rep("0",q), sd=rep("1",q) )
  # }
  # 
  # if(q>1){
  # V.raneff <- paste0("c(", paste0(raneff, collapse = ","), ")")
  # ran.loglike <- paste0("-0.5*",V.raneff, "%*%invSIGMA%*%", V.raneff,"+0.5*log(det(invSIGMA))-0.5*", q,"*log(2*pi)")
  # }
  #####
  
 
    if(is.null(nlmeObject$lower.disp))  nlmeObject$lower.disp <- rep(0, q)
    if(is.null(nlmeObject$upper.disp))  nlmeObject$upper.disp <- rep(Inf, q)
    
    lower.disp <- c(nlmeObject$lower.disp, sigma.lower.disp)
    upper.disp <- c(nlmeObject$upper.disp, sigma.upper.disp)
    names(lower.disp) <- names(upper.disp) <- c(disp,sigma.disp)
    
    if(is.null(nlmeObject$lower.fixed)) nlmeObject$lower.fixed <- rep(-Inf, p)
    if(is.null(nlmeObject$upper.fixed)) nlmeObject$upper.fixed <- rep(Inf, p)
    
    lower.fixed <- c(nlmeObject$lower.fixed, sigma.lower)
    upper.fixed <- c(nlmeObject$upper.fixed, sigma.upper)
    
    names(lower.fixed) <- names(upper.fixed) <- c(fixed, sigma.fixed)
  
    disp <- c(disp,sigma.disp)

    fixed <- c(fixed, sigma.fixed)

  
  if(!is.null(ranCovObject)){
    loglike <- list(mu.loglike=mu.loglike, sigma.loglike=sigma.loglike, 
                    randisp.loglike=randisp.loglike)
    return(list(loglike=loglike,
                fixed.par=fixed, ran.eff=c(raneff,sigma.raneff, randisp), disp.par=disp,
                str.fixed=str.fixed, str.disp=str.disp,
                lower.fixed=lower.fixed, upper.fixed=upper.fixed,
                lower.disp=lower.disp, upper.disp=upper.disp,
                response=resp,
                rvX=rvX,
                SIGMA.dim=q,
                ran.eff.dim=c(length(raneff),length(sigma.raneff), length(randisp)),
                sigma.df=sigma.df, randisp.df=ranCovObject$df,
                nf=nlmeObject$nf,
                muExpr=mu))
  } else {
    loglike <- list(mu.loglike=mu.loglike, sigma.loglike=sigma.loglike)
    return(list(loglike=loglike,
                fixed.par=fixed, ran.eff=c(raneff,sigma.raneff), disp.par=disp,
                str.fixed=str.fixed, str.disp=str.disp,
                lower.fixed=lower.fixed, upper.fixed=upper.fixed,
                lower.disp=lower.disp, upper.disp=upper.disp,
                response=resp,
                rvX=rvX,
                SIGMA.dim=q,
                ran.eff.dim=c(length(raneff),length(sigma.raneff)),
                sigma.df=sigma.df,
                nf=nlmeObject$nf,
                muExpr=mu))
  }
  
  
}
