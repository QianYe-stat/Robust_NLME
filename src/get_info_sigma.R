
get_info_sigma<- function(sigmaObejct){
  
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
    fixed <- paste("alpha", 0:(p-1), sep="") # name fixed pars
  } else {fixed=NULL}
  
  
  if(q>0){
    
    raneff <- paste("a", 0:(q-1), sep="")  # name random effects
    
    if(sigmaObject$ran.dist=="inverse-Chi"){
      
      raneff_loglike <- make_loglike_invChi(raneff, sigmaObject$df)
      
    }
  } else { 
    raneff <- NULL
    raneff_loglike <- NULL
  }
  
  linear_pred <- paste0(c(fixed, raneff), rep("*", p+q), 
                        c(rvX, rvZ), sep="", collapse="+")
  
  if(sigmaObject$link=="log"){
    sigma_expr <- paste0("exp(", linear_pred, ")")
  }
  
  str.val <- sigmaObject$str.val
  names(str.val) <- fixed
  
  return(list(sigmaExpr=sigma_expr, loglike=raneff_loglike, raneff=raneff, fixed=fixed, str.val=str.val))
}



