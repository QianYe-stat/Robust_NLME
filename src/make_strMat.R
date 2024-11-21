# This function returns a matrix with elements in string.
# Moreover, it returns L(l)'L(l)=SIGMA, which is a spherical
# parameterization with diagonal elements 1

make_strMat <- function(q2, block=NULL){
  
  if(is.null(block)){
    Mat <- make_Mat(q2)
    M <- Mat$M
    Mpar <- Mat$Mpar
  } else {
    k <- length(block)
    M.list <- vector("list",k)
    Mpar <- c()
    for(kk in 1:k){
      Mat <- make_Mat(block[kk], letters[kk])
      M.list[[kk]] <- Mat$M
      Mpar <- c(Mpar, Mat$Mpar)
    }
    
    cum_sum <- c(0, cumsum(block))
    num_block <- vector("list",k)
    
    for(kk in 1:k){
      num_block[[kk]] <- c((1+cum_sum[kk]):cum_sum[kk+1])
    }
    M <- matrix("0",q2,q2)  
    for(i in 1:q2){
      for(j in 1:q2){
        within <-c()
        for(kk in 1:k){
          within <- c(within, (i %in% num_block[[kk]]) & (j %in% num_block[[kk]]))
        }
        if(any(within)) {
          loc <- which(within)
          M[i,j] <- M.list[[loc]][i-cum_sum[loc], j-cum_sum[loc]]
        }
      }
    }
     
  }
 
  M.expr <- paste0("matrix(c(", paste0(c(M), collapse=","), "),nrow=", q2, ",ncol=", q2, 
                   ", byrow=TRUE)") 
  
  
 # dM <- array(NA, dim=c(q2,q2,Mp))
  
 # for(i in 1:Mp){
#    ttz <- unlist(lapply(M, function(x){Deriv(x, Mpar[i])}))
#    dM[,,i] <- matrix(as.character(ttz),q2,q2)
#  }
  
  return(list(M=M.expr, Mpar=Mpar))

# return(list(M=M.expr, Mpar=Mpar, dM=dM))
  
}