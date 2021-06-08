#' @param nlmeObject
#' @param long.data
#' @param idVar


long.data=dat
idVar="patid"

Rnlme <- function(nlmeObject, long.data, idVar, 
                  itertol=1e-3, Ptol=1e-2, iterMax=50, Verbose=FALSE)
##################################### settings for nlme model 
nlmeReturn <- get_nlme_loglike(nlmeObject)

Jloglike <- nlmeReturn$loglike # log likelihood conditional on random effects
Jraneff <- nlmeReturn$ran.eff # random effects 
Jfixed  <-  nlmeReturn$fixed.par   # fixed parameters 
Jdisp   <-  nlmeReturn$disp.par  # dispersion parameters
nf <- nlmeReturn$nf

p <- length(Jfixed)  #  dimension of fixed parameters 
q <- length(Jraneff) # dimension of random effects
q1 <- nlmeReturn$SIGMA.dim # dimension of SIGMA

##################################### initial values
fixedest0 <- nlmeReturn$str.fixed
dispest0 <- nlmeReturn$str.disp
invSIGMA0 <- SIGMA <- diag(1,q1,q1)
Lval0 <- NULL # re-parameters for COV matrix

#################################### bounds
lower.fixed <- NULL
lower.disp <- nlmeReturn$lower.disp
upper.disp <- nlmeReturn$upper.disp


#################################### settings for dataset
group <- long.data[ , idVar]  # grouping variable, e.g patient ID
uniqueID <- unique(group)   
n <- length(uniqueID)  # sample size
ni <- table(group)   # number of repeated measurements for each subject 
N <- nrow(long.data)


#################################### condition for iteration
likDiff <- Diff <- 1
convergence <- 1
M <- 1

itertol=1e-3
Ptol=1e-2
iterMax=50

while(likDiff > itertol & Diff > Ptol & M < iterMax) {
  #################################### estimation
  
  # estimate random effects
  ran.output <- est_raneff(RespLog=Jloglike, long.data, idVar, Jraneff, 
                           fixedest0, dispest0, invSIGMA0,
                           uniqueID, n,ni,q,N, Verbose=Verbose, scale=TRUE)
  Bi <- ran.output$Bi
  B <- ran.output$B
  print("estimate random effects --- done.")
  
  # estimate fixed parameters
  fixed.output <- est_fixed(RespLog=Jloglike, long.data, Jfixed,
                            fixedest0, dispest0, invSIGMA0,
                            Bi, B, 
                            lower=lower.fixed,
                            Verbose=Verbose)
  
  fixedest <- fixed.output$beta
  print("estimate fixed parameters --- done.")
  
  # estimate dispersion parameters
  disp.output <- est_dispersion(RespLog=Jloglike, long.data, Jdisp,
                                fixedest, dispest0, invSIGMA0, Lval0,
                                Bi, B,
                                lower=lower.disp, upper=upper.disp,
                                Verbose=Verbose)
  
  dispest <- disp.output$disp
  invSIGMA <- disp.output$invSIGMA
  Lval <- disp.output$Lval
  
  print("estimate dispersion parameters --- done.")
  
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
  
  # calcuate relative changes in mean parameters
  Diff <- mean(c(abs((fixedest - fixedest0)/(fixedest0 + 1e-6))))
  
  
  ############## print
  cat("############## Iteration:", M, "###############","\n")
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
  warning("Iteration limit reached without covergence.")
  convergence <- 1
}
if(likDiff <= itertol & likDiff >= 0){
  message("Successful convergence. Iteration stops because likDiff <= itertol.")
  convergence <- 0
}
if(Diff <= Ptol){
  message("Successful convergence. Iteration stops because FixedParDiff <= Ptol.")
  convergence <- 0
}

 
  # estimate sd's of parameter estimates  
  sd_output <- get_sd(RespLog=Jloglike, long.data,  
               fixedest0, dispest0, invSIGMA0,
               Bi, B,
               Jfixed, Jraneff)

  print("estimate SD for fixed parameters --- done.")
  
list(fixedest = fixedest0,
     fixedSD = sd_output,
     Bi = Bi, 
     #B = B,
     SIGMA = solve(invSIGMA0), 
     dispersion = dispest0,
     convergence = convergence==0,
     loglike_value = loglike_value
     #long.data = long.data,
     #surv.data = surv.data,
     #RespLog = RespLog,
     #idVar = idVar, uniqueID = uniqueID,
     #Jraneff = Jraneff
)
