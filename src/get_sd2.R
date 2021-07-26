get_sd2 <- function(RespLog, long.data, idVar,
                    fixedest0, dispest0, invSIGMA0,SIGMA0,
                    Bi, B, q_split,
                    Jfixed,Jraneff){
  
  n <- nrow(Bi)
  p <- length(Jfixed)
  
  
  dhat <- dispest0[c(1,2,3)]
  Dhat <- diag(dhat)%*%SIGMA0%*%diag(dhat)
  Tmat <- diag(0, p, p)
  
  for(i in 1:n){
    subdat <- subset(long.data, long.data[, idVar]==uniqueID[i])
    
    qL <- nrow(subdat)
    
    X <- Z <- as.matrix(cbind(rep(1, qL),subdat$day, subdat$cd4))
    Di <- Z %*% Dhat %*% t(Z)+diag(dispest0[4]^2, qL,qL)
    
    Tmat <- Tmat + t(X) %*% solve(Di) %*% X
    
  }
  sd <- sqrt(diag(solve(Tmat)))
  return(sd)
}