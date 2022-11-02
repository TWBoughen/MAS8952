brownian_kernel = function(x, y, amplitude = 1, ...) {
  # brownian motion kernel
  kern = function(x, y) {
    return(amplitude * pmin(x, y))
  }
  return(outer(x, y, kern))
}

squared_exponential_kernel = function(x, y, amplitude = 1, lengthscale = 1, ...) {
  # squared exponential kernel
  kern = function(x, y) {
    d = abs(x - y) / lengthscale
    return(amplitude * exp(- (d ** 2) / 2))
  }
  return(outer(x, y, kern))
}

exponential_kernel = function(x, y, amplitude = 1, lengthscale = 1, ...) {
  # exponential kernel
  kern = function(x, y) {
    d = abs(x - y) / lengthscale
    return(amplitude * exp(-d))
  }
  return(outer(x, y, kern))
}

matern_32_kernel = function(x, y, amplitude = 1, lengthscale = 1, ...) {
  # Matern kernel with nu = 3/2
  kern = function(x, y) {
    d = abs(x - y) / lengthscale
    term1 = 1 + sqrt(3) * d
    term2 = exp(-sqrt(3) * d)
    return(amplitude * term1 * term2)
  }
  return(outer(x, y, kern))
}

matern_52_kernel = function(x, y, amplitude = 1, lengthscale = 1, ...) {
  # Matern kernel with nu = 5/2
  kern = function(x, y) {
    d = abs(x - y) / lengthscale
    term1 = 1 + sqrt(5) * d + 5 / 3 * d ** 2
    term2 = exp(-sqrt(5) * d)
    return(amplitude * term1 * term2)
  }
  return(outer(x, y, kern))
}

rational_quadratic_kernel = function(x, y, amplitude = 1, lengthscale = 1, alpha = 1, ...) {
  # rational quadratic kernel
  kern = function(x, y) {
    d = abs(x - y) / lengthscale
    return(amplitude * ((1 + d ** 2 / (2 * alpha)) ** (-alpha)))
  }
  return(outer(x, y, kern))
}

periodic_kernel = function(x, y, amplitude = 1, lengthscale = 1, period = 1, ...) {
  # squared exponential kernel
  kern = function(x, y) {
    d = sin((pi / period) * abs(x - y)) / lengthscale
    return(amplitude * exp(- 2 * (d ** 2)))
  }
  return(outer(x, y, kern))
}

linear_kernel = function(x, y, amplitude = 1, c = 0, lambda_0 = 0, ...) {
  # linear kernel
  kern = function(x, y) {
    return(amplitude * ((x - c) * (y - c)) + lambda_0)
  }
  return(outer(x, y, kern))
}