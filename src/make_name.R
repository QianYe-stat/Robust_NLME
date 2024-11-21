make_name <- function(name, value){
  dat <- data.frame(as.list(value))
  names(dat) <- name
  return(dat)
}

