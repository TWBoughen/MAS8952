zero_mean = function(x, ...) {
  # zero mean function
  defaultW <- getOption("warn") 
  options(warn = -1)
  out=numeric(length(x))
  options(warn = defaultW)
  return(out)
  
  
}

constant_mean = function(x, c=1, ...) {
  # constant kernel
  defaultW <- getOption("warn") 
  options(warn = -1)
  out=numeric(length(x)) + c
  options(warn = defaultW)
  return(out)
}

linear_mean = function(x, a=1, b=1, ...) {
  # linear mean function
  return(a * x + b)
}