 est_raneff <- function(RespLog, long.data, idVar, Jraneff,
                       fixedest0, dispest0, invSIGMA0,
                       uniqueID, n,ni,q,N,q_split,df.sigma, df.randisp,
                       Verbose=TRUE, scale=TRUE){
  
  nBi <- nB <- c()
  
  for(i in 1:n){
    subdat <- subset(long.data, long.data[, idVar]==uniqueID[i])
    
    bi <- est_individual_raneff1(RespLog=RespLog, data=subdat, raneff=Jraneff, 
                                fixedest=fixedest0, dispest=dispest0, invSIGMA=invSIGMA0,
                                Verbose=Verbose)
    nBi <- rbind(nBi, bi)
    nB <- rbind(nB, matrix(rep(bi, ni[i]), ncol=q, byrow=T))           
    if(Verbose==TRUE)  cat("i=",i, "out of", n, "individuals.\n")    
  }
  if(scale==TRUE){
    
    # ### scale B1
    Bi <- as.matrix(nBi)%*%diag(1/apply(nBi,2,sd),q,q)
    B <- as.matrix(nB)%*%diag(1/apply(nBi,2,sd),q,q)
    cenBi <- apply(Bi, 2, mean)
    Bi <- Bi - matrix(rep(1, n),ncol=1)%*%matrix(cenBi, nrow=1)
    B <- B - matrix(rep(1, N),ncol=1)%*%matrix(cenBi, nrow=1)
    
    # q1 <- q_split[1]; q2 <- q_split[2]; q3 <- q_split[3]
    # 
    # nBi1 <- as.matrix(nBi[,c(1:q1)]); nB1 <- as.matrix(nB[,c(1:q1)])
    # nBi2 <- as.matrix(nBi[,-c(1:q1)][,c(1:q2)]); nB2 <- as.matrix(nB[,-c(1:q1)][,c(1:q2)])
    # nBi3 <- as.matrix(nBi[,-c(1:(q1+q2))]); nB3 <- as.matrix(nB[,-c(1:(q1+q2))])
    # 
    # ### scale B1
    # Bi1 <- as.matrix(nBi1)%*%diag(1/apply(nBi1,2,sd),q1,q1)
    # B1 <- as.matrix(nB1)%*%diag(1/apply(nBi1,2,sd),q1,q1)
    # cenBi1 <- apply(Bi1, 2, mean)
    # Bi1 <- Bi1 - matrix(rep(1, n),ncol=1)%*%matrix(cenBi1, nrow=1)
    # B1 <- B1 - matrix(rep(1, N),ncol=1)%*%matrix(cenBi1, nrow=1)
    
    ### scale B2
     
    # tg_mean <- log(df.sigma/(df.sigma-2))
    # tg_sd <- sqrt(2*df.sigma^2/(df.sigma-2)^2/(df.sigma-4))/(df.sigma/(df.sigma-2))
    # 
    # Bi2 <- as.matrix(nBi2)%*%diag(tg_sd/apply(nBi2,2,sd),q2,q2)
    # B2 <- as.matrix(nB2)%*%diag(tg_sd/apply(nBi2,2,sd),q2,q2)
    # cenBi2 <- apply(Bi2, 2, mean)
    # Bi2 <- Bi2 - matrix(rep(1, n),ncol=1)%*%matrix(cenBi2-tg_mean, nrow=1)
    # B2 <- B2 - matrix(rep(1, N),ncol=1)%*%matrix(cenBi2-tg_mean, nrow=1)
    # 
    # 
    # rm(tg_mean); rm(tg_sd)
    
    ### scale B3
    # tg_mean <- log(df.randisp/(df.randisp-2))
    # tg_sd <- sqrt(2*df.randisp^2/(df.randisp-2)^2/(df.randisp-4))/(df.randisp/(df.randisp-2))
    # 
    # Bi3 <- as.matrix(nBi3)%*%diag(tg_sd/apply(nBi3,2,sd),q3,q3)
    # B3 <- as.matrix(nB3)%*%diag(tg_sd/apply(nBi3,2,sd),q3,q3)
    # cenBi3 <- apply(Bi3, 2, mean)
    # Bi3 <- Bi3 - matrix(rep(1, n),ncol=1)%*%matrix(cenBi3-tg_mean, nrow=1)
    # B3 <- B3 - matrix(rep(1, N),ncol=1)%*%matrix(cenBi3-tg_mean, nrow=1)
  
    
    # Bi <- cbind(Bi1, nBi2, nBi3)
    # B <- cbind(B1, nB2, nB3)
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