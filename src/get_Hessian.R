
# calculate hessian matrix 

get_Hessian <- function(loglik, pars){
  loglik <- parse(text=loglik)
  q <- length(pars)
  result <- as.list(rep(NA, q^2))
  
  
  for(i in 1:q){
    for(j in 1:q){
      k <- (i-1)*q+j
      Fst <- Deriv(loglik, pars[i])
      if(suppressWarnings(str_detect(Fst, "as.matrix"))) {
        Fst <- suppressWarnings(str_remove(Fst, "as.matrix"))
        Fst <- parse(text=Fst)
      }
      result[[k]] <- Deriv(Fst, pars[j]) 
      names(result)[k] <- paste(pars[i], pars[j],sep=",")
    }
  }
  
  return(result)
}  
