#' @title Plot Age-specific Failure Age
#' @description This is age-specific failure rate plot.
#' @param x An object of class "nsr"
#' @param n The index of nest. For example, n=1 is the 1st nest.
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
#'  plot_q(x=out, 1)

plot_q <- function(x, n, xlab = "nest age", ylab = "individual age-specific nest failure rate", ...){
  plot(1:x$jj, x$q[n, 1:x$jj], type = "l", xlab = xlab, ylab = ylab)
}
