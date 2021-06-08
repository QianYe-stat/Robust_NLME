get_nlme_loglike <- function(nlmeObject){
  
  ############## Return from get_info_sigma #####################
  sigmaInfo <- get_info_sigma(nlmeObject$sigma)
  sigma.raneff <- sigmaInfo$raneff  # random effect in residual dispersion (residual random effects)
  sigma.fixed <- sigmaInfo$fixed
  sigma.loglike <- sigmaInfo$loglike
  sigma.str <- sigmaInfo$str.val
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
  
  sp <- strsplit(as.character(ranCovObject$varying.disp), "~",  fixed=T)[[2]]
  randisp.comp <- strsplit(sp, "+",  fixed=T)[[1]]
  randisp.comp <- str_trim(randisp.comp) 
  
  
  p <- length(fix.comp)  # dimension of fixed pars
  q <- length(ran.comp)  # dimension of random eff
  
  ran.ind <- c(1:p) * fix.comp %in% ran.comp
  randisp.ind <- c(1:p) * fix.comp %in% randisp.comp  # identify the random effects with varying dispersion
  
  if(p>0){
    fixed <- paste0("beta", 1:p) # name fixed pars
  } else {fixed=NULL}
  
  if(q>0){
    raneff <- c()
    disp <- c()
    randisp <- c()
    
    for(i in 1:p){
      sub <- ran.ind[i]
      if(sub!=0){
        raneff <- c(raneff,paste0("u", sub)) # name random effects
        disp <- c(disp, paste0("d", sub))    # name the fixed dispersion of random eff
      } else {
        raneff <- c(raneff, 0)
        disp <- c(disp, 0)
      }
      
      rm(sub)
      sub <- randisp.ind[i]
      if(sub!=0){
        randisp <- c(randisp, paste0("b", sub)) # name the random effects for varying dispersion of random effect
      } else {                                  # double random effects
        randisp <- c(randisp, 0)
      }
      
    }
  } else {
    raneff=NULL
    disp=NULL
    randisp=NULL}
  
  disp.eff <- paste0(disp, rep("*", p), paste0("exp(", randisp, ")"))
  toteff <- paste0(fixed, rep("+", p), paste0(raneff, rep("*",p), disp.eff)) 
  
  mu <- paste0("nf(", paste0(toteff, collapse = ","), ",day)")
  
  if (nlmeObject$family=="normal"){ 
    
    sigma <- sigmaInfo$sigmaExpr
    loglike <- paste0("-0.5*(", resp, "-", mu, 
                      ")^2/",sigma, "-0.5*log(",sigma,")-0.5*log(2*pi)")
    
  }
  
  raneff <- str_subset(raneff, "[^0]")
  disp <- str_subset(disp, "[^0]")
  randisp <- str_subset(randisp, "[^0]")
  
  if(ranCovObject$ran.dist=="inverse-Chi"){
    randisp.loglike <- make_loglike_invChi(randisp, ranCovObject$df)
  }
  
  #####
  str.fixed <- nlmeObject$str.fixed
  names(str.fixed) <- fixed
  
  str.disp <- c(nlmeObject$str.disp, sigma.str)
  names(str.disp) <- c(disp,sigma.fixed)
  
  V.raneff <- paste0("c(", paste0(raneff, collapse = ","), ")")
  ran.loglike <- paste0("-0.5*",V.raneff, "%*%invSIGMA%*%", V.raneff,"+0.5*log(det(invSIGMA))-0.5*", q,"*log(2*pi)")

  #####
  if(is.null(nlmeObject$lower.disp))  nlmeObject$lower.disp <- rep(0, q)
  if(is.null(nlmeObject$upper.disp))  nlmeObject$upper.disp <- rep(Inf, q)
  
  lower.disp <- c(nlmeObject$lower.disp, -Inf, -Inf)
  upper.disp <- c(nlmeObject$upper.disp, Inf, Inf)
  names(lower.disp) <- names(upper.disp) <- c(disp,sigma.fixed)
  
  return(list(loglike=list(mu.loglike=loglike, sigma.loglike=sigma.loglike, 
                           randisp.loglike=randisp.loglike,ran.loglike=ran.loglike),
              fixed.par=fixed, ran.eff=c(raneff,sigma.raneff, randisp), disp.par=c(disp,sigma.fixed),
              str.fixed=str.fixed, str.disp=str.disp,
              lower.disp=lower.disp, upper.disp=upper.disp,
              response=resp,
              rvX=rvX,
              SIGMA.dim=q,
              nf=nlmeObject$nf))
}