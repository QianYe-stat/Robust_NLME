est_dispersion <- function(RespLog, long.data, Jdisp,
                           fixedest, dispest0, invSIGMA0, Lval0,
                           Bi, B,
                           lower, upper,
                           Verbose=TRUE){
  
  n <- nrow(Bi)
  q1 <- length(Jdisp)
  q2 <- ncol(invSIGMA0)
  Lmat  <- make_strMat(q2) 
  qL <- length(Lmat$Mpar)
  ran.loglike <- str_replace_all(RespLog$ran.loglike, "invSIGMA", paste0("solve(",Lmat$M, ")"))
  
  Yrandisp <- !is.null(RespLog$randisp.loglike)
  Ysigma <- !is.null(RespLog$sigma.loglike)
  
  Jdisp_new <- c(Jdisp, Lmat$Mpar)
  
  # ff() returns the negative value of h-likelihood to be optimized.
  ff <- function(xx){  # xx are input values of dispersion parameters
    fy <- numeric(1)
    
    # assign values to the parameters
    par.val <- make_name(Jdisp_new , xx)  
    par.val <- c(fixedest, par.val)
    
    
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
    ran.val <- vector("list", n)
    for(i in 1:n){
      ran.val[[i]] <-  with(par.val, with(Bi[i,], eval(parse(text=ran.loglike))))
    }
    ran.val <- unlist(ran.val) 
    
    # evaluate the adjusted profile h-likelihood
    fy <- sum(mu.val)+sum(sigma.val)+sum(randisp.val)+sum(ran.val)
    
    return(-fy)
  }
  
  k <- length(Jdisp_new)
  
  gr.mu <- Deriv(RespLog$mu.loglike, Jdisp)
  
  gr.ran <- Deriv(ran.loglike, Jdisp_new)
  
  ############ gradient function
  
  gr <- function(xx){
    fy <- numeric(k)
    
    # assign values to the parameters
    par.val <- make_name(Jdisp_new , xx)  
    par.val <- c(fixedest, par.val)
    
    
    # derivative of sd parameters 
    gr.mu.val <- c()
    for(i in 1:nrow(B)){
      val <-  with(par.val, with(long.data[i,], with(B[i,], eval(parse(text=gr.mu)))))
      gr.mu.val <- rbind(gr.mu.val, val)
    }
    gr.mu.val <- c(apply(gr.mu.val, 2, sum), rep(0, qL))
    
    gr.ran.val <- c()
    for(i in 1:n){
      temp.val <-  with(par.val, with(Bi[i,], eval(parse(text=gr.ran))))
      gr.ran.val <- rbind(gr.ran.val, temp.val)
    }
    gr.ran.val <- as.vector(apply(gr.ran.val, 2, sum))
    
    fy <- gr.mu.val+gr.ran.val
    return(-fy)
  }
  
  str_val00 <- dispest0
  convge <- -1
  M <- 0
  
  if(Verbose==TRUE) check=1  else check=0
  
  # start iteration
  while(convge != 0 & M<50){
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





