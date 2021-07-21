

get_sd <- function(RespLog, long.data,
                   fixedest0, dispest0, invSIGMA0,
                   Bi, B, q_split,
                   Jfixed, Jraneff){
  p <- length(Jfixed)
  q <- length(Jraneff)
  q1 <- q_split[1]
  n <- nrow(Bi)
  
  pars <- c(Jfixed, Jraneff[c(1:q1)])

  
  Hmat.mu <- get_Hessian(RespLog$mu.loglike, pars)
 # Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, pars)
 # Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, pars)
  Hmat.ran <- get_Hessian(RespLog$ran.loglike, pars)
  
  par.val <- as.list(c(fixedest0, dispest0))
  par.val$invSIGMA <- invSIGMA0
  
  
  mu.val <- get_Hvalue(Hmat.mu, p+q1, long.data, par.val, B)
 # sigma.val <- get_Hvalue(Hmat.sigma, p+q, data=NULL, par.val, Bi)
 # randisp.val <- get_Hvalue(Hmat.randisp, p+q, data=NULL, par.val, Bi)
  ran.val <- bdiag(diag(0, p), -invSIGMA0*(n))
  
  #Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val+ran.val))
  Hval <-  as.matrix(-(mu.val+ran.val))
  covMat <- ginv(Hval)
  sd2 <- diag(covMat)[1:p]
  sd <- sqrt(sd2)
  return(sd)
}