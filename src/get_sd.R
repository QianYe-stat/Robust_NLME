get_sd <- function(RespLog, long.data,
                   fixedest0, dispest0, invSIGMA0,
                   Bi, B, q_split,
                   Jfixed, Jraneff){
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  
  p <- length(Jfixed)
  q <- length(Jraneff)
  q1 <- q_split[1]
  n <- nrow(Bi)
  
  pars <- c(Jfixed, Jraneff)
  
  
  Hmat.mu <- get_Hessian(RespLog$mu.loglike, pars)
  Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, pars)
  if(Yrandisp) Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, pars)
  # Hmat.ran <- get_Hessian(RespLog$ran.loglike, pars)
  
  par.val <- as.list(c(fixedest0, dispest0))
  par.val$invSIGMA <- invSIGMA0
  
  
  mu.val <- get_Hvalue(Hmat.mu, p+q, long.data, par.val, B)
  sigma.val <- get_Hvalue(Hmat.sigma, p+q, data=NULL, par.val, Bi)
  if(Yrandisp) randisp.val <- get_Hvalue(Hmat.randisp, p+q, data=NULL, par.val, Bi)
  ran.val <- bdiag(diag(0, p), -invSIGMA0*(n), diag(0, (q-q1)))
  
  
  
  if(Yrandisp) {
    Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val+ran.val))
  } else Hval <-  as.matrix(-(mu.val+sigma.val+ran.val))
  #Hval <-  as.matrix(-(mu.val+ran.val))
  #Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val))
  covMat <- ginv(Hval)
  sd2 <- diag(covMat)[1:p]
  sd <- sqrt(sd2)
  return(sd)
}