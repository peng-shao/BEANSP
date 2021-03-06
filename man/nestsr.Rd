% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nestsr.R
\name{nestsr}
\alias{nestsr}
\title{Nest Survival Rate Estimate}
\usage{
nestsr(jj, nx, nn, zl, zr, x, y, a, b, sigma, day, enc, covar, n0, ntotal,
  mi = FALSE, tha = FALSE, mutha = FALSE, vartha = FALSE, ul = FALSE,
  ur = FALSE)
}
\arguments{
\item{jj}{The number of time units(days) a nest is required to survive to be considered successful.}

\item{nx}{The number of covariates.}

\item{nn}{The number of observed nests (sample size).}

\item{zl}{The smallest possible number of time units from the first encounter date to the outcome date.}

\item{zr}{The largest possible number of time units from the first encounter date to the outcome date.}

\item{x}{The nest-specific covariate-matrix.}

\item{y}{The nest fate.}

\item{a}{The specified value for hyperparameter a. The prior gamma(a,b) is for age effect variances.}

\item{b}{The specified value for hyperparameter b. The prior gamma(a,b) is for age effect variances.}

\item{sigma}{The specified values for hyperparameters age-effect variances.}

\item{day}{The iniital values for age effect of outcome rates.}

\item{enc}{The initial values for age effect of encounter rates.}

\item{covar}{The initial values for the coefficients for covariates.}

\item{n0}{The number of burn-in cycles.}

\item{ntotal}{The number of total Gibbs cycles.}

\item{mi}{The indicator of missing(censored) nest fate. If fate is not missing,
the default value is FALSE}

\item{tha}{The missing probability coefficients. If fate is not missing,
the default value is FALSE}

\item{mutha}{The mean number of missing probability coefficients. If fate is not missing,
the default value is FALSE}

\item{vartha}{The variance number of missing probability coefficients. If fate is not missing,
the default value is FALSE}

\item{ul}{The youngest possible age that the nest could have been when first encounterd.
If fate is missing, the default value is FALSE.}

\item{ur}{The oldest possible age that the nest could have been when first encounterd.
If fate is missing, the default value is FALSE.}
}
\value{
An object of class "nsr" containing the following components:
\describe{
\item{\code{n0}}{The number of burn-in cycles.}
\item{\code{ntotal}}{The number of total Gibbs cycles.}
\item{\code{jj}}{nest period time.}
\item{\code{enc}}{numerical values of estimate of encounter age effect for all age.}
\item{\code{day}}{numerical values of estimate of outcome age effect for all age .}
\item{\code{sigma}}{numerical values of estimate of age effect variances.}
\item{\code{covar}}{numerical values of estimate of regression coefficients.}
\item{\code{q}}{numerical values of estimate of age-specific outcome rates.}
\item{\code{del}}{numerical values of estimate of age-specific encounter rates.}
\item{\code{sr}}{numerical values of estimate of individual age-specific survival rates.}
\item{\code{asr}}{numerical values of estimate of average age-specific survival rates.}
\item{\code{casr}}{numerical values of estimate of average cumulative age-specific survival rates.}
\item{\code{DIC}}{model selection criterion DIC value}
\item{\code{Dbar}}{a expectation measure of how well the model fits the data.}
\item{\code{pd}}{a measure for the effective number of parameters of the model.}
\item{\code{trace1}}{numerical trace values for encounter age effect for all age.}
\item{\code{trace2}}{numerical trace values for outcome age effect for all age.}
\item{\code{trace3}}{numerical trace values for regression coefficients.}
\item{\code{trace4}}{numerical trace values for age effect variances.}
\item{\code{venc}}{standard deviation of encounter age effect for all age.}
\item{\code{vday}}{standard deviation of outcome age effect for all age.}
\item{\code{vsigma}}{standard deviation of age effect variances.}
\item{\code{vcovar}}{standard deviation of regression coefficients.}
\item{\code{vq}}{standard deviation of age-specific outcome rates.}
\item{\code{vdel}}{standard deviation of age-specific outcome rates.}
\item{\code{vsr}}{standard deviation of individual age-specific survival rates.}
\item{\code{vasr}}{standard deviation of average age-specific survival rates.}
\item{\code{vcasr}}{standard deviation of average cumulative age-specific survival rates.}
\item{\code{tha}}{numerical values of estimate of \code{theta}.}
\item{\code{mutha}}{numerical values of expectation of \code{theta}.}
\item{\code{varthat}}{standard deviation of \code{theta}.}
\item{\code{trace5}}{numerical trace values for \code{theta}}
\item{\code{WAIC}}{model selection criterion WAIC value}
}
}
\description{
This function can handle nests with unknown nest age, nest-specific
covariates(both discrete and continues covariates), and any hazard rate function(or
equivelently survival rate function) as long as it is a smooth function.
}
\details{
The Bayesian estimate of parameter is computed from its posterior
distribution which is simulated by Gibbs sampler. Users need to specify a set
of initial values ,the number of burn-in cycles and the total number of Gibbs
sampling cycles.
}
\examples{

n0 <- as.integer(10)
ntotal <- as.integer(110)
a <- as.double(rep(2.0,2))
b <- as.double(rep(1.0,2))
sigma <- as.double(rep(1.0,2))

data("example1")
jj <- as.integer(19)
nx <- as.integer(2)
nn <- as.integer(217)
x1 <- example1[1:nn,7]
x2 <- example1[1:nn,8]
x <- cbind(as.double(x1),as.double(x2))
temp <- nestconv(0,nn,jj,example1, censored  =  FALSE)
ul <- as.integer(temp[,2])
ur <- as.integer(temp[,3])
zl <- as.integer(temp[,4])
zr <- as.integer(temp[,5])
y <- as.integer(temp[,6])

day <- as.double(rep(0.0,jj))
enc <- as.double(rep(0.0,jj - 1))
covar <- as.double(rep(0.0,nx))
out <- nestsr(jj = jj, nx = nx, nn = nn, zl = zl, zr = zr, x = x, y = y,
              a = a, b = b, sigma = sigma, day = day, enc = enc, covar = covar,
              n0 = n0, ntotal = ntotal, ul = ul, ur = ur)

data("example4")
jj <- as.integer(26)
nx <- as.integer(1)
nn <- as.integer(221)
x <- as.double(as.numeric(example4[,7]) - 1)
temp <- nestconv(1,nn,jj,example4, censored  =  TRUE)
zl <- as.integer(temp[,2])
zr <- as.integer(temp[,3])
y <- as.integer(temp[,4])
mi <- as.integer(example4[,6])

day <- as.double(rep(0.0,jj))
enc <- as.double(rep(0.0,jj - 1))
covar <- as.double(rep(0.0,nx))
tha <- as.double(rep(0.0,2))
mutha <- as.double(rep(0.0,2))
vartha <- as.double(rep(5.0,2))
out <- nestsr(jj = jj, nx = nx, nn = nn, zl = zl, zr = zr, x = x, y = y, a = a,
              b = b, sigma = sigma, day = day, enc = enc, covar = covar, n0 = n0,
              ntotal = ntotal, mi = mi, tha = tha, mutha = mutha, vartha = vartha)

}
\references{
Cao, J., He, C., Suedkamp Wells, K.M., Millspaugh, J.J., and
Ryan, M.R. (2009).  Modeling age and nest-specific survival using a hierarchical
Bayesian approach. Biometrics, 65, 1052-1062.
}
\author{
Chong He, Yiqun Yang, Jing Cao, Peng Shao
}
\keyword{methods}
