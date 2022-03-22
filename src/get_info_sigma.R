
get_info_sigma<- function(sigmaObject){
  
  sp1 <- strsplit(as.character(as.formula(sigmaObject$model)), "~",  fixed=T)
  
  if( "(" %in% strsplit(as.character(sp1[[2]]), "")[[1]] ){  # if model includes random effects 
    sps <- strsplit(as.character(sp1[[2]]), "(",  fixed=T)[[1]]
    sp2 <- as.character(sps[1])
    rvX <- strsplit(sp2,"+",fixed=T)[[1]] 
    rvX <- rvX[-length(rvX)]   # returns covariates of fixed parameters
    
    sp3 <- as.character(strsplit(sps[2], ")",  fixed=T)[[1]])
    # returns covariates of random effects
    rvZ <- strsplit(as.character(strsplit(sp3,"|",fixed=T)[[1]][1]), "+", fixed=T)[[1]]           
    
    
  } else {
    rvX <- strsplit(as.character(sp1[[3]]), "+", fixed=T)[[1]]
    rvZ <- NULL
  }
  
  p <- length(rvX)  # dimension of fixed pars
  q <- length(rvZ) # dimension of random effects
  
  if(p>0){
    fixed <- paste(sigmaObject$fixName, 0:(p-1), sep="") # name fixed pars
  } else {fixed=NULL}
  
  
  if(q>0){
    
    raneff <- paste(sigmaObject$ranName, 0:(q-1), sep="")  # name random effects

    
    if(sigmaObject$ran.dist=="inverse-Chi"){
      df <- sigmaObject$df
      
      raneff_loglike <- make_loglike_invChi(raneff, df)
      
      linear_pred <- paste0(c(fixed, raneff), rep("*", p+q), 
                            c(rvX, rvZ), sep="", collapse="+")
      
    } else if(sigmaObject$ran.dist=="normal"){
      
      disp.par <- paste(sigmaObject$dispName, 0:(q-1), sep="")
      
      raneff_loglike <- make_loglike_normal(raneff, mean=rep("0",q), sd=rep("1",q) )
      
      linear_pred <- paste0(c(fixed, paste0(raneff, rep("*",q), disp.par)), rep("*", p+q),
                            c(rvX, rvZ), sep="", collapse="+")
      df<- NULL
    } else if(sigmaObject$ran.dist=="stdnormal"){
      
      raneff_loglike <- make_loglike_normal(raneff, mean=rep("0",q), sd=rep("1",q) ) 
      
      linear_pred <- paste0(c(fixed, raneff), rep("*", p+q), 
                            c(rvX, rvZ), sep="", collapse="+")
      df<- NULL
      
    }
  } else { 
    raneff <- NULL
    raneff_loglike <- NULL
  }
  

  if(sigmaObject$link=="log"){
    sigma_expr <- paste0("exp(", linear_pred, ")")
  }
  
  
  str.fixed <- sigmaObject$str.fixed
  lower.fixed <- sigmaObject$lower.fixed
  upper.fixed <- sigmaObject$upper.fixed
  

  if(is.null(lower.fixed)) lower.fixed <- rep(-Inf, p)
  if(is.null(upper.fixed)) upper.fixed <- rep(Inf, p)
  
  names(str.fixed) <- names(lower.fixed) <- names(upper.fixed) <- fixed
  disp <- str.disp <- lower.disp <- upper.disp <- NULL
  
  if(sigmaObject$ran.dist=="normal"){
    
  str.disp <- sigmaObject$str.disp
  lower.disp <- sigmaObject$lower.disp
  upper.disp <- sigmaObject$upper.disp
  
  if(is.null(lower.disp)) lower.disp <- rep(0, q)
  if(is.null(upper.disp)) upper.disp <- rep(Inf, q)
  
  names(str.disp) <- names(lower.disp) <- names(upper.disp) <- disp.par
  
  disp <- disp.par
  
  }
  
  trueVal.model <- sigmaObject$trueVal.model
  if(!is.null(trueVal.model)){
    model<- trueVal.model$model
    true.expr <- get_mu(model)
    
    sigma_expr <- str_replace_all(sigma_expr, trueVal.model$var , paste0("(", true.expr , ")"))
  }
  
  return(list(sigmaExpr=sigma_expr, loglike=raneff_loglike, raneff=raneff, fixed=fixed, 
              df=df,str.val=str.fixed,lower=lower.fixed, upper=upper.fixed, 
              disp=disp, str.disp=str.disp, lower.disp=lower.disp, upper.disp=upper.disp ))
}



