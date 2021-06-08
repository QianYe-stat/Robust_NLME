est_raneff <- function(RespLog, long.data, idVar, Jraneff,
                       fixedest0, dispest0, invSIGMA0,
                       uniqueID, n,ni,q,N, Verbose=TRUE, scale=TRUE){
  
  nBi <- nB <- c()
  
  for(i in 1:n){
    subdat <- subset(long.data, long.data[, idVar]==uniqueID[i])
    
    bi <- est_individual_raneff(RespLog=RespLog, data=subdat, raneff=Jraneff, 
                                fixedest=fixedest0, dispest=dispest0, invSIGMA=invSIGMA0,
                                Verbose=Verbose)
    nBi <- rbind(nBi, bi)
    nB <- rbind(nB, matrix(rep(bi, ni[i]), ncol=q, byrow=T))           
    if(Verbose==TRUE)  cat("i=",i, "out of", n, "individuals.\n")    
  }
  if(scale==TRUE){
    Bi <- as.matrix(nBi)%*%diag(1/apply(nBi,2,sd),q,q)
    B <- as.matrix(nB)%*%diag(1/apply(nBi,2,sd),q,q)
    cenBi <- apply(Bi, 2, mean)
    Bi <- Bi - matrix(rep(1, n),ncol=1)%*%matrix(cenBi, nrow=1)
    B <- B - matrix(rep(1, N),ncol=1)%*%matrix(cenBi, nrow=1)
  } else {
    Bi <- nBi
    B <- nB
  }
  
  Bi <- as.data.frame(Bi)
  B <- as.data.frame(B)
  names(B) <- names(Bi) <- Jraneff
  rownames(Bi) <- rownames(B) <- c()
  
  return(output=list(Bi=Bi, B=B))
}