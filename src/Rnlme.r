#' @param nlmeObject
#' @param long.data
#' @param idVar


Rnlme <- function(nlmeObject, long.data, idVar, 
                  itertol=1e-3, Ptol=1e-2, iterMax=20, Verbose=FALSE){
  #set.seed(123)
  ##################################### settings for nlme model 
  nlmeReturn <- get_nlme_loglike(nlmeObject)
  
  Jloglike <- nlmeReturn$loglike # log likelihood conditional on random effects
  Jraneff <- nlmeReturn$ran.eff # random effects 
  Jfixed  <-  nlmeReturn$fixed.par   # fixed parameters 
  Jdisp   <-  nlmeReturn$disp.par  # dispersion parameters
  nf <- nlmeReturn$nf
  
  p <- length(Jfixed)  #  dimension of fixed parameters 
  q <- length(Jraneff) # dimension of random effects
  q_split <-  nlmeReturn$ran.eff.dim
  qSIGMA <- nlmeReturn$SIGMA.dim # dimension of SIGMA
  
  df.sigma <- nlmeReturn$sigma.df
  df.randisp <- nlmeReturn$randisp.df
  
  ##################################### initial values
  fixedest0 <- nlmeReturn$str.fixed
  dispest0 <- nlmeReturn$str.disp
  invSIGMA0 <- SIGMA <- diag(1,qSIGMA,qSIGMA)
  Lval0 <- NULL # re-parameters for COV matrix
  
  #################################### bounds
  lower.fixed <- nlmeReturn$lower.fixed
  upper.fixed <- nlmeReturn$upper.fixed
  lower.disp <- nlmeReturn$lower.disp
  upper.disp <- nlmeReturn$upper.disp
  
  
  #################################### settings for dataset
  group <- long.data[ , idVar]  # grouping variable, e.g patient ID
  uniqueID <- unique(group)   
  n <- length(uniqueID)  # sample size
  ni <- table(group)   # number of repeated measurements for each subject 
  N <- nrow(long.data)
  
  
  #################################### condition for iteration
  likDiff <- Diff <- Diff0 <- 1
  convergence <- 1
  M <- 1
  
  while(!(likDiff <= itertol | (Diff <= Ptol & Diff0 <= Ptol) | M > iterMax)) {
    #################################### estimation
    Diff0 <- Diff
    cat("############## Iteration:", M, "###############","\n")
    
    # estimate random effects
    cat("Start estimating random effects ...\n")
    ran.output <- est_raneff(RespLog=Jloglike, long.data, idVar, Jraneff, 
                             fixedest0, dispest0, invSIGMA0,
                             uniqueID, n,ni,q, N, q_split,df.sigma, df.randisp,
                             Verbose=Verbose, scale=TRUE)
    Bi <- ran.output$Bi
    B <- ran.output$B
    cat("done.\n")
    
    # estimate fixed parameters
    cat("Start estimating fixed parameters ... \n")
    fixed.output <- est_fixed(RespLog=Jloglike, long.data, Jfixed,
                              fixedest0, dispest0, invSIGMA0,
                              Bi, B, 
                              lower=lower.fixed, upper=upper.fixed,
                              Verbose=Verbose)
    
    fixedest <- fixed.output$beta
    cat("done.\n")
    
    # estimate dispersion parameters
    cat("Start estimating dispersion parameters ... \n")
    disp.output <- est_dispersion(RespLog=Jloglike, long.data, Jdisp,
                                  fixedest, dispest0, invSIGMA0, Lval0,
                                  Bi, B,
                                  lower=lower.disp, upper=upper.disp,
                                  Verbose=Verbose)
    
    dispest <- disp.output$disp
    invSIGMA <- disp.output$invSIGMA
    Lval <- disp.output$Lval
    
    cat("done.\n")
    
    ####################################################    
    ################## update results ##################
    ####################################################
    # calculate approximated log marginal likelihood value
    loglike_value <- get_loglike_value(RespLog=Jloglike, long.data,
                                       fixedest, dispest, invSIGMA, Bi, B)
    
    
    if(M == 1){
      likDiff <- 1
    } else{
      likDiff <- abs(loglike_value-loglike_value0)/abs(loglike_value0)
    }
    
    # calculate relative changes in mean parameters
    Diff <- mean(c(abs((fixedest - fixedest0)/(fixedest0 + 1e-6))))
    
    
    ############## print
    
    cat("fixed.par:", round(fixedest, 2), "\n")
    cat("FixedParDiff = ", Diff, '\n')
    cat("likDiff = ", likDiff, '\n')
    if(!is.null(unlist(dispest))){
      cat("dispersion.par:", round(unlist(dispest), 2), "\n")
    }
    cat("loglike:", loglike_value, "\n")
    cat("##########################################","\n")
    
    fixedest0 <- fixedest
    dispest0 <- dispest
    invSIGMA0 <- invSIGMA
    Lval0 <- Lval
    loglike_value0 <- loglike_value
    M <- M+1  
  }
  
  ## messages about convergence success or failure
  if((likDiff > itertol  & Diff>Ptol ) ){
    message("Iteration limit reached without covergence.")
    convergence <- 1
  }
  if(likDiff <= itertol & likDiff >= 0){
    message("Successful convergence. Iteration stops because likDiff <= itertol.")
    convergence <- 0
  }
  if(Diff <= Ptol & Diff0 <= Ptol){
    message("Successful convergence. Iteration stops because FixedParDiff <= Ptol in two consective iteration.")
    convergence <- 0
  }
  
  cat("Start estimating SD for fixed parameters ...\n ...\n")
  
  # estimate sd's of parameter estimates  
  sd_output <- get_sd(RespLog=Jloglike, long.data,  
                      fixedest0, dispest0, invSIGMA0,
                      Bi, B, q_split,
                      Jfixed, Jraneff)
  
  cat("done.\n")
  
  #### AIC
  AIC <- 2*length(c(fixedest0, dispest0, Lval0)) -2*loglike_value0
  BIC <- length(c(fixedest0, dispest0, Lval0))*log(nrow(B))-2*loglike_value0
  
  list(fixedest = fixedest0,
       fixedSD = sd_output,
       Bi = Bi, 
       #B = B,
       SIGMA = solve(invSIGMA0), 
       dispersion = dispest0,
       convergence = convergence==0,
       loglike_value = loglike_value0,
       AIC=AIC,
       BIC=BIC
       #long.data = long.data,
       #surv.data = surv.data,
       #RespLog = RespLog,
       #idVar = idVar, uniqueID = uniqueID,
       #Jraneff = Jraneff
  )
}
