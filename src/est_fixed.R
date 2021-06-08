est_fixed <- function(RespLog, long.data, Jfixed,
                      fixedest0, dispest0, invSIGMA0,
                      Bi, B, 
                      lower,
                      Verbose=TRUE){
  
  n <- nrow(Bi)
  
  # ff() returns the negative value of the h-likelihood to be optimized.
  
  ff <- function(xx){
    fy <- numeric(1)
    
    # assign values to parameters
    par.val <- make_name(Jfixed, xx)    
    par.val <- c(par.val, dispest0)
    par.val$invSIGMA <- invSIGMA0
    
    # evaluate h-likelihood
    mu.val <- with(long.data, with(par.val, with(B, eval(parse(text=RespLog$mu.loglike)))))
    sigma.val <-  with(par.val, with(Bi, eval(parse(text=RespLog$sigma.loglike))))
    randisp.val <- with(par.val, with(Bi, eval(parse(text=RespLog$randisp.loglike))))
    
    ran.val <- vector("list", n)
    for(i in 1:n){
      ran.val[[i]] <-  with(par.val, with(Bi[i,], eval(parse(text=RespLog$ran.loglike))))
    }
    ran.val <- unlist(ran.val)
    
    # evaluate profile h-likelihood
    fy <- sum(mu.val)+sum(sigma.val)+sum(randisp.val)+sum(ran.val)  
    
    return(-fy)
  }
  
  k <-  length(Jfixed)
  
  gr.mu <- Deriv(RespLog$mu.loglike, Jfixed)
  
  
  gr <- function(xx){
    fy <- numeric(k)
    
    # assign values to parameters
    par.val <- make_name(Jfixed, xx)    
    par.val <- c(par.val, dispest0)
    par.val$invSIGMA <- invSIGMA0
    
    
    gr.mu.val <- gr.sigma.val <- gr.randisp.val <- gr.ran.val <- rep(NA, k)
    
    ##
    val <- with(long.data, with(par.val, with(B, eval(parse(text=gr.mu)))))
    val <- matrix(val, ncol=k, byrow=FALSE)
    gr.mu.val <- as.vector(apply(val, 2, sum))
    
    fy <-  gr.mu.val
    
    return(-fy)
  }
  
  ## start iteration
  str_val0 <- fixedest0
  convge <- -1
  M <- 0
  
  if(Verbose==TRUE) check<-1  else check<-0
  
  if(is.null(lower)){
    lower= -Inf; method="BFGS"
  } else {
    method="L-BFGS-B"
  }
  
  while(convge != 0 & M<20){
    
    str_val0 <- sapply(str_val0, function(x)x+rnorm(1,0, min(1, abs(x/5))))
    
    result <- try(optim(par=str_val0, fn=ff, gr=gr,method=method,lower=lower,
                        control = list(trace=check,maxit=1000)),silent=TRUE)
    
    error_mess <- attr(result, "class")
    
    if(is.null(error_mess)){ 
      convge <- result$convergence 
      str_val0 <- result$par
    } else {
      convge = -1
    }
    M <- M+1
    if(Verbose==TRUE){ 
      cat(paste0("\n M=",M,", Convergence=", convge==0,".\n"))
      print(result)}
    
  }
  
  if(convge==0){
    beta <- result$par
    names(beta) <- Jfixed
    fval <- result$value
  } else {
    stop("Iteration stops because fixed parameters can not be successfully estimated.")  
  }
  
  return(list(beta=beta, fval=fval))
  
}

