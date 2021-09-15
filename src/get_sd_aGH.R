get_sd_aGH <- function(RespLog, long.data, idVar, 
                       fixedest0, dispest0, invSIGMA0,Bi, B,
                       Jfixed, Jraneff,  
                       ghsize=4, Silent=T, epsilon=10^{-6}, 
                       parallel=F){
  q <- ncol(Bi)
  group <- long.data[ , idVar]  
  uniqueID <- unique(group)   
  
  GHzsamp0 = mgauss.hermite(n=ghsize, mu=rep(0,q), sigma=NULL)
  
  idSIGMA = get_idSIGMA_aGH(RespLog, long.data, idVar, uniqueID,
                            fixedest0, dispest0, invSIGMA0,Bi, B,
                            Jfixed, Jraneff) 
  # generate GH samples by subject
  GHsample0 <- as.list(rep(NA,n))
  for(i in 1:n){
    GHsample0[[i]] = mgauss.hermite(n=ghsize, mu=as.numeric(Bi[i,]), sigma=idSIGMA[[i]])
  }
  
  GHsd2 = try(calculate_aGH(RespLog, long.data, idVar, uniqueID,
                            fixedest0, dispest0, invSIGMA0,
                            GHzsamp0,GHsample0,
                            Jfixed, Jraneff,
                            ghsize, epsilon, parallel),  silent = Silent)
  return(GHsd2)
}