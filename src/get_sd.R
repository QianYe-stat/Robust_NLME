

get_sd <- function(RespLog, long.data,
                   fixedest0, dispest0, invSIGMA0,
                   Bi, B,
                   Jfixed, Jraneff){
  p <- length(Jfixed)
  q <- length(Jraneff)
  
  pars <- c(Jfixed, Jraneff)
  
  Hmat.mu <- get_Hessian(RespLog$mu.loglike, pars)
  Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, pars)
  Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, pars)
  Hmat.ran <- get_Hessian(RespLog$ran.loglike, pars)
  
  par.val <- as.list(c(fixedest0, dispest0))
  par.val$invSIGMA <- invSIGMA0
  
  
  mu.val <- get_Hvalue(Hmat.mu, p+q, long.data, par.val, B)
  sigma.val <- get_Hvalue(Hmat.sigma, p+q, data=NULL, par.val, Bi)
  randisp.val <- get_Hvalue(Hmat.randisp, p+q, data=NULL, par.val, Bi)
  ran.val <- get_Hvalue(Hmat.ran, p+q, data=NULL, par.val, Bi)
  
  Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val+ran.val))
  covMat <- ginv(Hval)
  sd2 <- diag(covMat)[1:p]
  sd <- sqrt(sd2)
  return(sd)
}