get_sd1 <- function(RespLog, long.data, idVar,
                    fixedest0, dispest0, invSIGMA0,SIGMA0,
                    Bi, B, q_split,
                    Jfixed, Jraneff){
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  p <- length(Jfixed)
  q <- length(Jraneff)
  q1 <- q_split[1]
  n <- nrow(Bi)
  
  pars <- c(Jfixed, Jraneff[c(1:q1)])
  
  
  Hmat.mu <- get_Hessian(RespLog$mu.loglike, pars)
  if(Ysigma) Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, pars)
  if(Yrandisp) Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, pars)
  # Hmat.ran <- get_Hessian(RespLog$ran.loglike, pars)
  
  par.val <- as.list(c(fixedest0, dispest0))
  par.val$invSIGMA <- invSIGMA0
  
  
  mu.val <- get_Hvalue(Hmat.mu, p+q1, long.data, par.val, B)
  
  if(Ysigma){
    sigma.val <- get_Hvalue(Hmat.sigma, p+q1, data=NULL, par.val, Bi)
  } else sigma.val <- diag(0, p+q1, p+q1)
  
  if(Yrandisp) {
    randisp.val <- get_Hvalue(Hmat.randisp, p+q1, data=NULL, par.val, Bi)
  } else randisp.val <- diag(0, p+q1, p+q1)
  
  ran.val <- bdiag(diag(0, p,p), -invSIGMA0*(n))
  
  
  Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val+ran.val))
  
  #Hval <-  as.matrix(-(mu.val+ran.val))
  #Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val))
  covMat <- solve(Hval)
  sd2 <- diag(covMat)[1:p]
  sd <- sqrt(sd2)
  return(sd)
}