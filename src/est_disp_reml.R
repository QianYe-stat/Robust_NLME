est_disp_reml <- function(RespLog, long.data, Jdisp, Jfixed, Jraneff,
                           fixedest, dispest0, invSIGMA0, Lval0,
                           Bi, B,
                           lower, upper,
                           Verbose=TRUE){
  
  n <- nrow(Bi)
  q1 <- length(Jdisp)
  p <- length(Jfixed)
  q <- length(Jraneff)
  q2 <- ncol(invSIGMA0)
  Lmat  <- make_strMat(q2) 
  qL <- length(Lmat$Mpar)
  Mpar <- Lmat$Mpar
  Bi_sub <- Bi[,c(1:q2)]
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  Jdisp_new <- c(Jdisp, Mpar)
  
  
  ## Hessian matrix
  Dpars <- c(Jfixed, Jraneff)
  
  qD <- length(Dpars)
  
  Hmat.mu <- get_Hessian(RespLog$mu.loglike, Dpars)
  if(Ysigma) Hmat.sigma <- get_Hessian(RespLog$sigma.loglike, Dpars)
  if(Yrandisp) Hmat.randisp <- get_Hessian(RespLog$randisp.loglike, Dpars)
  
  
  # ff() returns the negative value of h-likelihood to be optimized.
  ff <- function(xx){  # xx are input values of dispersion parameters
    fy <- numeric(1)
    
    # assign values to the parameters
    par.val <- make_name(Jdisp,xx[1:q1])  
    par.val <- c(fixedest, par.val)
    
    Lval <- make_name(Lmat$Mpar, xx[-(1:q1)])
    
    SIGMA <- with(Lval, eval(parse(text=Lmat$M)))
    invSIGMA <- solve(SIGMA)
    
    
    # evaluate the h-likelihood
    mu.val <- with(long.data, with(par.val, with(B, eval(parse(text=RespLog$mu.loglike)))))
    if(Ysigma){
      sigma.val <- with(par.val, with(Bi, eval(parse(text=RespLog$sigma.loglike))))
    } else sigma.val <- 0
    
    if(Yrandisp) {
      randisp.val <- with(par.val, with(Bi, eval(parse(text=RespLog$randisp.loglike))))
    } else {
      randisp.val <- 0
    }
    
    ran.val <- 0.5*log(det(invSIGMA)) - diag(0.5*as.matrix(Bi_sub)%*%invSIGMA%*%t(as.matrix(Bi_sub)))
    
    
    ########### Evaluate the H matrix #################
    
    Hmu.val <- get_Hvalue(Hmat.mu, qD, long.data, par.val, B)
    
    
    if(Ysigma){
      Hsigma.val <- get_Hvalue(Hmat.sigma, qD, data=NULL, par.val, Bi)
    } else Hsigma.val <- diag(0, qD, qD)
    
    if(Yrandisp) {
      Hrandisp.val <- get_Hvalue(Hmat.randisp, qD, data=NULL, par.val, Bi)
    } else Hrandisp.val <- diag(0, qD, qD)
    
    
    Hran.val <-  bdiag(diag(0,p), -invSIGMA*n, diag(0, q-q2))
    
    Hval <-  as.matrix(-(Hmu.val+Hsigma.val+Hrandisp.val+Hran.val))
    
    # evaluate the adjusted profile h-likelihood
    fy <- sum(mu.val)+sum(sigma.val)+sum(randisp.val)+sum(ran.val)-0.5*log(det(Hval/2/pi))
    
    return(-fy)
    
  }
  
  k <- length(Jdisp_new)
  
  gr.mu <- Deriv(RespLog$mu.loglike, Jdisp)
  
  dH.mu <- as.list(rep(NA, qD^2))
  
  for(i in 1:qD^2){
    dH.mu[[i]] <- Deriv(Hmat.mu[[i]] , Jdisp)
  }
  
  
  ############ gradient function
  
  gr <- function(xx){
    fy <- numeric(k)
    
    # assign values to the parameters
    par.val <- make_name(Jdisp,xx[1:q1])  
    par.val <- c(fixedest, par.val)
    
    Lval <- make_name(Lmat$Mpar, xx[-(1:q1)])
    
    SIGMA <- with(Lval, eval(parse(text=Lmat$M)))
    invSIGMA <- solve(SIGMA)
    
    ########### Evaluate the H matrix #################
    Hmu.val <- get_Hvalue(Hmat.mu, qD, long.data, par.val, B)
    
    if(Ysigma){
      Hsigma.val <- get_Hvalue(Hmat.sigma, qD, data=NULL, par.val, Bi)
    } else Hsigma.val <- diag(0, qD, qD)
    
    if(Yrandisp) {
      Hrandisp.val <- get_Hvalue(Hmat.randisp, qD, data=NULL, par.val, Bi)
    } else Hrandisp.val <- diag(0, qD, qD)
    
    Hran.val <-  bdiag(diag(0,p), -invSIGMA*n, diag(0, q-q2))
    
    Hval <-  as.matrix(-(Hmu.val+Hsigma.val+Hrandisp.val+Hran.val))
    
    invH <-  solve(Hval)
    
    ########## derivative of disp parameters ############### 
    gr.mu.val <- c()
    for(i in 1:nrow(B)){
      val <-  with(par.val, with(long.data[i,], with(B[i,], eval(parse(text=gr.mu)))))
      gr.mu.val <- rbind(gr.mu.val, val)
    }
    gr.mu.val <- apply(gr.mu.val, 2, sum)
    
    dH.mu.val <- matrix(NA, nrow=qD^2, ncol=q1)
    for(i in 1:qD^2){
      val <- c()
      for(j in 1:nrow(B)){
        val <- rbind(val,with(par.val, with(long.data[j,], with(B[j,], eval(parse(text=dH.mu[[i]]))))))
      }
      dH.mu.val[i,] <- apply(val, 2,sum)
      rm(val)
    }
    
    traces.mu <- rep(NA, q1)
    for(i in 1:q1){
      dH.mu.disp <- matrix(dH.mu.val[,i], qD,qD)
      traces.mu[i] <-  tr(invH%*%(-dH.mu.disp)) 
    }
    
    fy1 <- -(gr.mu.val-0.5*traces.mu)
    
    ################ derivative of disp_new parameters ###############
    
    fy2 <- rep(NA, qL)
    for(i in 1:qL){
      dM <- Lmat$dM[,,i]
      dM_val <- get_Hvalue(dM, q2, par.val=Lval)
      gh1 <- sum(-0.5* tr(invSIGMA%*%dM_val)+
                   0.5*diag(as.matrix(Bi_sub)%*%invSIGMA%*%dM_val%*%invSIGMA%*%t(as.matrix(Bi_sub))))
      DH <- invSIGMA%*%dM_val%*%invSIGMA*n
      DH2 <- bdiag(diag(0,p), DH, diag(0,(q-q2))) 
      fy2[i] <-  -(gh1-0.5*tr(as.matrix(invH%*%(-DH2))))
    }
    
    
    fy <- c(fy1, fy2)
    return(fy)
  }
  
  str_val00 <- dispest0
  convge <- -1
  M <- 0
  
  if(Verbose==TRUE) check=1  else check=0
  
  # start iteration
  while(convge != 0 & M<1){
    cat("\n", M)
    if(is.null(Lval0)){
      Lval00 <- runif( qL, 0, pi)
    } else{
      Lval00 <- Lval0+rnorm(qL,0,0.01)
    }
    str_val0 <-  unlist(c(str_val00+rnorm(q1,0,0.1), Lval00))
    
    
    result <- try(optim(par=str_val0, fn=ff, gr=gr,
                        method="L-BFGS-B",
                        lower=c(lower, rep(0.01,q2)), 
                        upper=c(upper,rep(pi*0.95,q2)),
                        control = list(trace=check,maxit=1000)), 
                  silent=T)
    
    error_mess <- attr(result, "class")
    
    if(is.null(error_mess)){ 
      convge <- result$convergence 
      str_val00 <- result$par[1:q1]
      Lval0 <- result$par[-c(1:q1)]
    } else {
      str_val00 <- dispest0
      convge = -1
      
    }
    M <- M+1
    if(Verbose==TRUE){ 
      cat(paste0("\n M=",M,", Convergence=", convge==0,".\n"))
      print(result)}
  }
  
  if(convge==0){
    if(q1>0){
      output <- result$par[1:q1]
      names(output) <- Jdisp
      Lval <- make_name(Lmat$Mpar, result$par[-(1:q1)])
      SIGMAest <- with(Lval, eval(parse(text=Lmat$M)))
      invSIGMAest <- solve(SIGMAest)
    } else{ 
      output=NULL
      Lval <- make_name(Lmat$Mpar, result$par)
      SIGMAest <- with(Lval, eval(parse(text=Lmat$M)))
      invSIGMAest <- solve(SIGMAest)
    }
    
  } else {
    stop("Iteration stops because dispersion parameters can not be successfully estimated.")  
  }
  
  return(list(disp=output, invSIGMA=invSIGMAest, Lval=Lval, SIGMA=SIGMAest))
}





