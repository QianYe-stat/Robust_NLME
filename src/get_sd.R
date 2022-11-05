get_sd <- function(RespLog, long.data, idVar,
                   fixedest0, dispest0, invSIGMA0,SIGMA0,
                   Bi, B,
                   Jfixed, Jraneff){
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  p <- length(Jfixed)
  q <- length(Jraneff)
  #q1 <- q_split[1]
  n <- nrow(Bi)

  pars <- c(Jfixed, Jraneff)
  
  Hmat.mu <- map(RespLog$mu.loglike, function(t){
     get_Hessian(t, pars)
  })
  
  if(Ysigma) Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, pars)
  if(Yrandisp) Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, pars)
  Hmat.ran <- get_Hessian(RespLog$ran.loglike, pars)

  
  par.val <- as.list(c(fixedest0, dispest0))
  par.val$invSIGMA <- invSIGMA0
  
  
  mu.val <- map(Hmat.mu, function(t){
    get_Hvalue(t, p+q, long.data, par.val, B)
  })
 
  mu.val <- Reduce('+', mu.val)
  
  if(Ysigma){
    sigma.val <- get_Hvalue(Hmat.sigma, p+q, data=NULL, par.val, Bi)
  } else sigma.val <- diag(0, p+q, p+q)
  
  if(Yrandisp) {
    randisp.val <- get_Hvalue(Hmat.randisp, p+q, data=NULL, par.val, Bi)
  } else randisp.val <- diag(0, p+q, p+q)
  
  #ran.val <- bdiag(diag(0, p,p), -invSIGMA0*(n), diag(0, (q-q1),(q-q1)))
  ran.val <- get_Hvalue(Hmat.ran, p+q, data=NULL, par.val, Bi)
  
  Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val+ran.val))
  
  #Hval <-  as.matrix(-(mu.val+ran.val))
  #Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val))
  covMat <- solve(Hval)

  sd2 <- diag(covMat)[1:p]
  sd <- sqrt(sd2)
  return(sd)
}