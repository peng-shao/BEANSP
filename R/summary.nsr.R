#' @title Summary \code{nsr} Object
#' @rdname summary.nsr
#' @aliases print.summary.nsr
#' @rdname summary.nsr
#' @description Convert date format data or interval format data to the standard format.
#' @param object An object of class "nsr".
#' @param covnames A character vector of the names of covariates.
#' @param x An object of class "summary.nsr".
#' @param ... additional arguments affecting the summary produced.
#' @return An object of class "summary.nsr" containing the following components:
#' \describe{
#' \item{\code{covar.table}}{A data frame of mean, standard deviation, 2.5% quantile
#' and 97.5 quantile of covariate(s).}
#' \item{\code{DIC}}{Model selection criterion DIC value}
#' \item{\code{WAIC}}{Model selection criterion WAIC value}
#' \item{\code{n0}}{The number of burn-in cycles.}
#' \item{\code{ntotal}}{The number of total Gibbs cycles.}
#' }
#' @author
#' Chong He, Yiqun Yang, Jing Cao, Peng Shao
#' @keywords methods
#' @importFrom stats quantile sd
#' @export
#' @examples
#'
#' data("example4")
#' n0 <- as.integer(10)
#' ntotal <- as.integer(110)
#' a <- as.double(rep(2.0,2))
#' b <- as.double(rep(1.0,2))
#' sigma <- as.double(rep(1.0,2))
#' jj <- as.integer(26)
#' nx <- as.integer(1)
#' nn <- as.integer(221)
#' x <- as.double(as.numeric(example4[,7]) - 1)
#' temp <- nestconv(1,nn,jj,example4, censored  =  TRUE)
#' zl <- as.integer(temp[,2])
#' zr <- as.integer(temp[,3])
#' y <- as.integer(temp[,4])
#' mi <- as.integer(example4[,6])
#'
#' day <- as.double(rep(0.0,jj))
#' enc <- as.double(rep(0.0,jj - 1))
#' covar <- as.double(rep(0.0,nx))
#' tha <- as.double(rep(0.0,2))
#' mutha <- as.double(rep(0.0,2))
#' vartha <- as.double(rep(5.0,2))
#' out <- nestsr(jj = jj, nx = nx, nn = nn, zl = zl, zr = zr, x = x, y = y, a = a,
#'               b = b, sigma = sigma, day = day, enc = enc, covar = covar, n0 = n0,
#'               ntotal = ntotal, mi = mi, tha = tha, mutha = mutha, vartha = vartha)
#' summary(out)
#'
#' data("example1")
#' jj <- as.integer(19)
#' nx <- as.integer(2)
#' nn <- as.integer(217)
#' x1 <- example1[1:nn,7]
#' x2 <- example1[1:nn,8]
#' x <- cbind(as.double(x1),as.double(x2))
#' temp <- nestconv(0,nn,jj,example1, censored  =  FALSE)
#' ul <- as.integer(temp[,2])
#' ur <- as.integer(temp[,3])
#' zl <- as.integer(temp[,4])
#' zr <- as.integer(temp[,5])
#' y <- as.integer(temp[,6])
#'
#' day <- as.double(rep(0.0,jj))
#' enc <- as.double(rep(0.0,jj - 1))
#' covar <- as.double(rep(0.0,nx))
#' out <- nestsr(jj = jj, nx = nx, nn = nn, zl = zl, zr = zr, x = x, y = y,
#'               a = a, b = b, sigma = sigma, day = day, enc = enc, covar = covar,
#'               n0 = n0, ntotal = ntotal, ul = ul, ur = ur)
#' summary(out)
#'

summary.nsr <- function(object, covnames = FALSE, ...){
  covtabl <- function(object){
    return(c(Estimate = mean(object),  SD = sd(object),
             `2.5% quantile` = unname(quantile(object, probs = 0.025)),
             `97.5% quantile` = unname(quantile(object, probs = 0.975))))
  }
  covar.trace <- as.data.frame(object$trace3)
  if (identical(covnames, FALSE)) {
    colnames(covar.trace) <- paste("Cov", 1:ncol(covar.trace))
  }else{
    if (ncol(covar.trace) != length(covnames)) {
      stop("number of covariate names does not equal number of columns of covariate trace.")
    }
    colnames(covar.trace) <- covnames
  }
  covar.table <- as.data.frame(t(apply(covar.trace, 2, covtabl)))
  DIC <- object$DIC
  WAIC <- object$WAIC
  ntotal <- object$ntotal
  n0 <- object$n0
  res <- list(covar.table = covar.table, DIC = DIC, WAIC = WAIC, n0 = n0, ntotal = ntotal)
  class(res) <- "summary.nsr"
  return(res)
}

#' @rdname summary.nsr
#' @export

print.summary.nsr <- function(x, ...){
  cat(c("Interations:", x$ntotal, ", Burn-in:", x$n0, "\n\n Estimation Table\n"))
  print(x$covar.table)
  cat(c("\nInformation Criterion:\n DIC:",x$DIC))
  if (!is.na(x$WAIC)){
    cat(c("\n WAIC:",x$WAIC))
  }
  cat("\n")
}

