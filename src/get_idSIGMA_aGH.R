

get_idSIGMA_aGH <- function(RespLog, long.data, idVar, uniqueID,
                           fixedest0, dispest0, invSIGMA0,Bi, B,
                           Jfixed, Jraneff){ 
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  q <- length(Jraneff)
  n <- nrow(Bi)
  
  Hmat.mu <- get_Hessian(RespLog$mu.loglike, pars=Jraneff)
  if(Ysigma) Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, pars=Jraneff)
  if(Yrandisp) Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, pars=Jraneff)
  Hmat.ran <- get_Hessian(RespLog$ran.loglike, pars=Jraneff)
  
  
  par.val <- as.list(c(fixedest0, dispest0))
  par.val$invSIGMA <- invSIGMA0
  
  idSIGMA = as.list(rep(NA,n))
  
  for(i in 1:n){
    indexi <- long.data[,idVar]==uniqueID[i]
    subdat <- subset(long.data, indexi)
   
    mu.val <- get_Hvalue(Hmat.mu, q, subdat, par.val, raneff.val=B[indexi,])
    
    if(Ysigma){
      sigma.val <- get_Hvalue(Hmat.sigma, q, data=NULL, par.val, Bi[i,])
    } else sigma.val <- diag(0, q, q)
    
    if(Yrandisp) {
      randisp.val <- get_Hvalue(Hmat.randisp, q, data=NULL, par.val, Bi[i,])
    } else randisp.val <- diag(0, q, q)
  
    ran.val <- get_Hvalue(Hmat.ran, q, data=NULL, par.val, raneff.val=Bi[i,])
    
    H <-  -(mu.val+sigma.val+randisp.val+ran.val)

    idSIGMA[[i]] <- solve(H)
  }
  
  return(idSIGMA)
  
}



