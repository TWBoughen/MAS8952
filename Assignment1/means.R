zero_mean = function(x, ...) {
  # zero mean function
  return(numeric(length(x)))
}

constant_mean = function(x, c=1, ...) {
  # constant kernel
  return(numeric(length(x)) + c)
}

linear_mean = function(x, a=1, b=1, ...) {
  # linear mean function
  return(a * x + b)
}
