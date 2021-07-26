get_sd1 <- function(RespLog, long.data,
                   fixedest0, dispest0, invSIGMA0,
                   Bi, B, q_split,
                   Jfixed, Jraneff){
  p <- length(Jfixed)
  q <- length(Jraneff)
  q1 <- q_split[1]
  n <- nrow(Bi)
  
  pars <- c(Jfixed)
  
  
  Hmat.mu <- get_Hessian(RespLog$mu.loglike, pars)
 
  par.val <- as.list(c(fixedest0, dispest0))
  par.val$invSIGMA <- invSIGMA0
  
  
  mu.val <- get_Hvalue(Hmat.mu, p, long.data, par.val, B, n=n)

  #ran.val <- bdiag(diag(0, p), -invSIGMA0*(n))
  
  #Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val+ran.val))
  Hval <-  as.matrix(-(mu.val))
  covMat <- ginv(Hval)
  sd2 <- diag(covMat)[1:p]
  sd <- sqrt(sd2)
  return(sd)
}