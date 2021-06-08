get_loglike_value <- function(RespLog, long.data, fixedest, dispest, invSIGMA, Bi, B){
  
  par.val <- as.list(c(fixedest, dispest))
  par.val$invSIGMA <- invSIGMA
  
  mu.val <- with(long.data, with(par.val, with(B,eval(parse(text=RespLog$mu.loglike)))))
  sigma.val <- with(par.val, with(Bi,eval(parse(text=RespLog$sigma.loglike))))
  randisp.val <- with(par.val, with(Bi,eval(parse(text=RespLog$randisp.loglike))))
  
  n <- nrow(Bi)
  ran.val <- vector("list", n)
  for(i in 1:n){
    ran.val[[i]] <-  with(par.val, with(Bi[i,], eval(parse(text=RespLog$ran.loglike))))
  }
  ran.val <- unlist(ran.val)
  
  return(sum(mu.val)+sum(sigma.val)+sum(randisp.val)+sum(ran.val))
}