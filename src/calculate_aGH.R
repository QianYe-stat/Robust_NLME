calculate_aGH <- function(RespLog, long.data, idVar, uniqueID,
                          fixedest0, dispest0, invSIGMA0,
                          GHzsamp0,GHsample0,
                          Jfixed, Jraneff,
                          ghsize, epsilon, parallel){
  
  

  #############
  
  p <- length(Jfixed)
  q <- length(Jraneff)
  
  gr.mu <- Deriv(RespLog$mu.loglike, Jfixed)
  
  weights <- GHzsamp0$weights
  n <- length(uniqueID)
  
 
  ## function for calculating S(theta, bi)
  gr <- function(xx){
    fy <- numeric(p)
    # assign values to parameters
    par.val <- make_name(Jfixed, xx)
    par.val <-  c(par.val, dispest0)
    par.val$invSIGMA <- invSIGMA0
   
    gn = matrix(NA, nrow=n, ncol=p)
    
    for(i in 1:n){
      # i=1
      subdat <-  subset(long.data, long.data[,idVar]==uniqueID[i])
    
      Bi_nodes <-  as.data.frame(GHsample0[[i]]$points)
      names(Bi_nodes) <- Jraneff
      
      likefn =  rep(NA, ghsize^q)
      gri = matrix(NA, ghsize^q, p)
      
      norm_term = 0
      
      for(j in 1:ghsize^q){
        # j=1
        
        llike.val <- map(RespLog, function(tt){
          with(subdat, with(par.val, with(Bi_nodes[j,], eval(parse(text=tt)))))
        })
        
        likeli_all <- exp(sum(unlist(llike.val)))
        
        likefn[j] = likeli_all*exp(sum(GHzsamp0$points[j,]^2))
        norm_term = norm_term+likeli_all
        
        gr.mul.val <- with(subdat, with(par.val, with(Bi_nodes[j,], eval(parse(text=gr.mu)))))
        
        gri[j,] = c(apply(matrix(gr.mul.val, byrow=FALSE, ncol=p), 2, sum))*likefn[j]
      }
      
      gn[i,] = as.vector(t(gri)%*%weights)/c(norm_term)
    }
    
    fy = apply(gn,2,sum)*det(invSIGMA0)^{-1/2}*2^{q/2}
    return(-fy)
  }
  
  ## esimate Hessian matrix by numerical derivative 
  est = unlist(c(fixedest0))
  
  if(parallel==TRUE){
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    
    Hmat = foreach(exponent = 1:(p), 
                   .combine = rbind
    )  %dopar%  
      {
        setwd(here::here("src"))
        file.sources = list.files(pattern="*.R$")
        sapply(file.sources,source,.GlobalEnv)
        
        Delta = rep(0, p)
        Delta[exponent] = epsilon
        gr1= gr(est+Delta)
        gr2= gr(est-Delta)
        (gr1-gr2)/(epsilon*2)
      }
    
    stopCluster(cl)
  } else {
    Hmat = matrix(NA, nrow=p, ncol=p)
    for(i in 1:(p)){
      # i=1
      Delta = rep(0, p)
      Delta[i] = epsilon
      Hmat[i,]=(gr(est+Delta)-gr(est))/epsilon
      # cat("i=",i,'\n')
    }
  }
  
  sd = sqrt(diag(solve(Hmat)))
  
  return(sd)
}