#' @title Trace Plot for Variance Effect
#' @description Displays a plot of sampled values for the variance effect
#' vs. iterations.
#' @param x An object of class "nsr"
#' @param j The j-th age effect.j=1,2,...,jj-1
#' @param n0 The number of burn-in cycles.
#' @param ntotal The number of total Gibbs cycles.
#' @param xlab Label of \code{x}-axis
#' @param ylab Label of \code{y}-axis
#' @param ... Arguments to be passed to methods.
#' @importFrom graphics plot
#' @author
#' Chong He, Yiqun Yang, Jing Cao, Peng Shao
#' @keywords methods
#' @export
#' @examples
#'
#'  n0 <- as.integer(10)
#'  ntotal <- as.integer(110)
#'  a <- as.double(rep(2.0,2))
#'  b <- as.double(rep(1.0,2))
#'  sigma <- as.double(rep(1.0,2))
#'
#'  data("example1")
#'  jj <- as.integer(19)
#'  nx <- as.integer(2)
#'  nn <- as.integer(217)
#'  x1 <- example1[1:nn,7]
#'  x2 <- example1[1:nn,8]
#'  x <- cbind(as.double(x1),as.double(x2))
#'  temp <- nestconv(0,nn,jj,example1, censored  =  FALSE)
#'  ul <- as.integer(temp[,2])
#'  ur <- as.integer(temp[,3])
#'  zl <- as.integer(temp[,4])
#'  zr <- as.integer(temp[,5])
#'  y <- as.integer(temp[,6])
#'
#'  day <- as.double(rep(0.0,jj))
#'  enc <- as.double(rep(0.0,jj - 1))
#'  covar <- as.double(rep(0.0,nx))
#'  out <- nestsr(jj = jj, nx = nx, nn = nn, zl = zl, zr = zr, x = x, y = y,
#'                a = a, b = b, sigma = sigma, day = day, enc = enc, covar = covar,
#'                n0 = n0, ntotal = ntotal, ul = ul, ur = ur)
#'  plot_var_trace(x=out, j=1, n0=n0, ntotal=ntotal)

plot_var_trace <- function(x, j, n0, ntotal, xlab = "iter", ylab = paste("encounter age effect at j=",j), ...){
  iter <- (n0 + 1):ntotal
  trace <- x$trace4[iter, j]
  plot(iter, trace, type = "l", xlab = xlab, ylab = ylab)
}
