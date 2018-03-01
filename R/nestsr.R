#' @title Nest Survival Rate Esitmate
#' @description This function can handle nests with unknown nest age, nest-specific
#' covariates(both discrete and continues covariates), and any hazard rate function(or
#' equivelently survival rate function) as long as it is a smooth function.
#' @param jj The number of time units(days) a nest is required to survive to be considered successful.
#' @param nx The number of covariates.
#' @param nn The number of observed nests (sample size).
#' @param zl The smallest possible number of time units from the first encounter date to the outcome date.
#' @param zr The largest possible number of time units from the first encounter date to the outcome date.
#' @param x The nest-specific covariate-matrix.
#' @param y The nest fate.
#' @param a The specified value for hyperparameter a. The prior gamma(a,b) is for age effect variances.
#' @param b The specified value for hyperparameter b. The prior gamma(a,b) is for age effect variances.
#' @param sigma The specified values for hyperparameters age-effect variances.
#' @param day The iniital values for age effect of outcome rates.
#' @param enc The initial values for age effect of encounter rates.
#' @param covar The initial values for the coefficients for covariates.
#' @param n0 The number of burn-in cycles.
#' @param ntotal The number of total Gibbs cycles.
#' @param mi The indicator of missing(censored) nest fate.
#' @param tha The missing probability coefficients.
#' @param mutha The mean number of missing probability coefficients.
#' @param vartha The variance number of missing probability coefficients.
#' @param ul The youngest possible age that the nest could have been when first encounterd.
#' @param ur The oldest possible age that the nest could have been when first encounterd.
#' @details The Bayesian estimate of parameter is computed from its posterior
#' distribution which is simulated by Gibbs sampler. Users need to specify a set
#' of initial values ,the number of burn-in cycles and the total number of Gibbs
#' sampling cycles.
#' @return An object of class "nsr" containing the following components:
#' \describe{
#' \item{\code{n0}}{The number of burn-in cycles.}
#' \item{\code{ntotal}}{The number of total Gibbs cycles.}
#' \item{\code{jj}}{nest period time.}
#' \item{\code{enc}}{numerical values of estimate of encounter age effect for all age.}
#' \item{\code{day}}{numerical values of estimate of outcome age effect for all age .}
#' \item{\code{sigma}}{numerical values of estimate of age effect variances.}
#' \item{\code{covar}}{numerical values of estimate of regression coefficients.}
#' \item{\code{q}}{numerical values of estimate of age-specific outcome rates.}
#' \item{\code{del}}{numerical values of estimate of age-specific encounter rates.}
#' \item{\code{sr}}{numerical values of estimate of individual age-specific survival rates.}
#' \item{\code{asr}}{numerical values of estimate of average age-specific survival rates.}
#' \item{\code{casr}}{numerical values of estimate of average cumulative age-specific survival rates.}
#' \item{\code{DIC}}{model selection criterion DIC value}
#' \item{\code{Dbar}}{a expectation measure of how well the model fits the data.}
#' \item{\code{pd}}{a measure for the effective number of parameters of the model.}
#' \item{\code{trace1}}{numerical trace values for encounter age effect for all age.}
#' \item{\code{trace2}}{numerical trace values for outcome age effect for all age.}
#' \item{\code{trace3}}{numerical trace values for regression coefficients.}
#' \item{\code{trace4}}{numerical trace values for age effect variances.}
#' \item{\code{venc}}{standard deviation of encounter age effect for all age.}
#' \item{\code{vday}}{standard deviation of outcome age effect for all age.}
#' \item{\code{vsigma}}{standard deviation of age effect variances.}
#' \item{\code{vcovar}}{standard deviation of regression coefficients.}
#' \item{\code{vq}}{standard deviation of age-specific outcome rates.}
#' \item{\code{vdel}}{standard deviation of age-specific outcome rates.}
#' \item{\code{vsr}}{standard deviation of individual age-specific survival rates.}
#' \item{\code{vasr}}{standard deviation of average age-specific survival rates.}
#' \item{\code{vcasr}}{standard deviation of average cumulative age-specific survival rates.}
#' }
#' @useDynLib BEANSP, .registration = TRUE
#' @importFrom stats var
#' @author
#' Chong He, Yiqun Yang, Jing Cao, Peng Shao
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
#'  out <- nestsr(jj = jj, nx = nx,nn = nn,zl = zl,zr = zr,x = x,y = y,a = a,b = b,
#'                sigma = sigma,day = day,enc = enc,covar = covar,n0 = n0,ntotal = ntotal,
#'                ul = ul, ur = ur)
#'
#'  data("example4")
#'  jj <- as.integer(26)
#'  nx <- as.integer(1)
#'  nn <- as.integer(221)
#'  x <- as.double(as.numeric(example4[,7]) - 1)
#'  temp <- nestconv(1,nn,jj,example4, censored  =  TRUE)
#'  zl <- as.integer(temp[,2])
#'  zr <- as.integer(temp[,3])
#'  y <- as.integer(temp[,4])
#'  mi <- as.integer(example4[,6])
#'
#'  day <- as.double(rep(0.0,jj))
#'  enc <- as.double(rep(0.0,jj - 1))
#'  covar <- as.double(rep(0.0,nx))
#'  tha <- as.double(rep(0.0,2))
#'  mutha <- as.double(rep(0.0,2))
#'  vartha <- as.double(rep(5.0,2))
#'  out <- nestsr(jj = jj,nx = nx,nn = nn,zl = zl,zr = zr,x = x,y = y,a = a,b = b,
#'                sigma = sigma, day = day,enc = enc,covar = covar,n0 = n0,ntotal = ntotal,
#'                mi = mi, tha = tha, mutha = mutha, vartha = vartha)
#'

nestsr <-  function(jj,nx,nn,zl,zr,x,y,a,b,sigma,day,enc,covar,n0,ntotal,
                    mi = FALSE, tha = FALSE, mutha = FALSE, vartha = FALSE,
                    ul = FALSE, ur = FALSE){
  if ((!identical(mi, FALSE)) && (!identical(tha, FALSE)) &&
      (!identical(mutha, FALSE)) && (!identical(vartha, FALSE))){
    out <- .Fortran("mbeansp",jj=jj,nx=nx,nn=nn,zl=zl,zr=zr,x=x,y=y,a=a,b=b,sigma=sigma,
                    day=day,enc=enc,tha=tha,covar=covar,n0=n0,ntotal=ntotal,
                    q=mat.or.vec(nn,jj+1),del=double(jj),DIC=double(1),Dbar=double(1),
                    pd=double(1),trace1=mat.or.vec(ntotal,jj-1),trace2=mat.or.vec(ntotal,jj),
                    trace3=mat.or.vec(ntotal,nx),trace4=mat.or.vec(ntotal,2),venc=double(jj-1),
                    vday=double(jj),vsigma=double(2),vcovar=double(nx),
                    vdel=double(jj),vq=mat.or.vec(nn,jj+1),asr=double(jj),vasr=double(jj),
                    casr=double(jj),vcasr=double(jj),sr1=mat.or.vec(nn,jj),vsr=mat.or.vec(nn,jj),
                    mutha=mutha,vartha=vartha,trace5=mat.or.vec(ntotal,2),
                    mi=mi,WAIC=double(1))
    vtha <- apply(out$trace5[(n0 + 1):ntotal,], 2, var)
    res <- list(n0=out$n0,ntotal=out$ntotal,jj=out$jj,zl=out$zl,zr=out$zr,y=out$y,
                enc=out$enc,day=out$day,sigma=out$sigma,covar=out$covar,q=out$q,
                del=out$del,WAIC=out$WAIC,DIC=out$DIC,Dbar=out$Dbar,pd=out$pd,trace1=out$trace1,
                trace2=out$trace2,trace3=out$trace3,trace4=out$trace4,venc=out$venc,
                vday=out$vday,vsigma=out$vsigma,vcovar=out$vcovar,vdel=out$vdel,
                vq=out$vq,asr=out$asr,vasr=out$vasr,casr=out$casr,vcasr=out$vcasr,
                sr=out$sr1,vsr=out$sr,tha=out$tha,mutha=out$mutha,vartha=out$vartha,
                trace5=out$trace5,ul=NA,ur=NA)
  }else{if ((!identical(ul, FALSE)) && (!identical(ur, FALSE))){
    out <- .Fortran("BEANSP",jj=jj,nx=nx,nn=nn,ul=ul,ur=ur,zl=zl,zr=zr,x=x,y=y,a=a,b=b,sigma=sigma,
                    day=day,enc=enc,covar=covar,n0=n0,ntotal=ntotal,q=mat.or.vec(nn,jj+1),
                    del=double(jj),DIC=double(1),Dbar=double(1),pd=double(1),
                    trace1=mat.or.vec(ntotal,jj-1),trace2=mat.or.vec(ntotal,jj),
                    trace3=mat.or.vec(ntotal,nx),trace4=mat.or.vec(ntotal,2),
                    venc=double(jj-1),vday=double(jj),vsigma=double(2),vcovar=double(nx),
                    vdel=double(jj),vq=mat.or.vec(nn,jj+1),asr=double(jj),vasr=double(jj),
                    casr=double(jj),vcasr=double(jj),sr=mat.or.vec(nn,jj),vsr=mat.or.vec(nn,jj))
    res <- list(n0=out$n0,ntotal=out$ntotal,jj=out$jj,ul=out$ul,ur=out$ur,zl=out$zl,
                zr=out$zr,y=out$y,enc=out$enc,day=out$day,sigma=out$sigma,covar=out$covar,
                q=out$q,del=out$del,DIC=out$DIC,Dbar=out$Dbar,pd=out$pd,trace1=out$trace1,
                trace2=out$trace2,trace3=out$trace3,trace4=out$trace4,venc=out$venc,
                vday=out$vday,vsigma=out$vsigma,vcovar=out$vcovar,vdel=out$vdel,vq=out$vq,
                asr=out$asr,vasr=out$vasr,casr=out$casr,vcasr=out$vcasr,sr=out$sr,vsr=out$sr,
                tha=NA,mutha=NA,vartha=NA,WAIC=NA,trace5=NA)
  }else{
    stop("Inputs are not consistent")
  }}
  class(res) <- "nsr"
  return(res)
}

