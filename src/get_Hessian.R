
# calculate hessian matrix 

get_Hessian <- function(loglik, pars){
  loglik_expr <- parse(text=loglik)
  q <- length(pars)
  result <- as.list(rep(NA, q^2))
  
  
  for(i in 1:q){
    for(j in 1:q){
      k <- (i-1)*q+j
      
      Fst_expr <- Deriv(loglik_expr, pars[i])
      Fst_char <- Deriv(loglik, pars[i])
      
      if(any(str_detect(Fst_char, "as.matrix"))) {
        
        Fst_char <- str_remove(Fst_char, "as.matrix")
        result[[k]] <- Deriv(Fst_char, pars[j])
      } else {
        result[[k]] <- Deriv(Fst_expr, pars[j]) 
      }
      names(result)[k] <- paste(pars[i], pars[j],sep=",")
    }
  }
  
  return(result)
}  
