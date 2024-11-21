get_sd_dipsersion <- function(RespLog, long.data, idVar,
                   fixedest0, dispest0, invSIGMA0,SIGMA0, Lval0,Lmat,
                   Bi, B,
                   Jfixed, Jraneff, Jdisp){
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  n <- nrow(Bi)
  q2 <- ncol(invSIGMA0)

 
  Lpar <- Lmat$Mpar
  ran.loglike <- str_replace_all(RespLog$ran.loglike, "invSIGMA", paste0("solve(",Lmat$M, ")"))
  
  pars <- c(Jdisp, Lpar)
  
  p <- length(pars)
  p1 <- length(Jdisp)
  p2 <- length(Lpar)
  

  Hmat.mu <- map(RespLog$mu.loglike, function(t){
    get_Hessian(t, Jdisp)
  })
  #if(Ysigma) Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, pars)
  #if(Yrandisp) Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, pars)
  if(!is.null(Lpar)) Hmat.ran <- get_Hessian(ran.loglike, Lpar)
  
  par.val <- as.list(c(fixedest0, dispest0, Lval0))
  
  mu.val <- map(Hmat.mu, function(t){
    get_Hvalue(t, p1, long.data, par.val, B)
  })
  mu.val <- Reduce('+', mu.val)
  
  mu.val <- bdiag(mu.val, diag(0, p2,p2))
  
  ran.val <- diag(0,p2,p2)
  
  if(!is.null(Lpar)){
    for(i in 1:n){
      temp.val <- get_Hvalue(Hmat.ran, p2, data=NULL, par.val, Bi[i,])
      ran.val <- ran.val+temp.val
    }
  }
  
  ran.val <- bdiag(diag(0, p1,p1), ran.val)
  
  
  Hval <-  as.matrix(-(mu.val+ran.val))
  
  #Hval <-  as.matrix(-(mu.val+ran.val))
  #Hval <-  as.matrix(-(mu.val+sigma.val+randisp.val))
  covMat <- solve(Hval)
  sd2 <- diag(covMat)[1:p1]
  sd <- sqrt(sd2)
  return(sd)
}