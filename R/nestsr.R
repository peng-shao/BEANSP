#' @title A test
#' @description a test function.
#' @param x this is a number
#' @return a
#' @return b
#' @useDynLib BEANSP, .registration = TRUE
#' @importFrom stats var
#' @export

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
    vtha <- apply(out$trace5[(n0+1):ntotal,], 2, var)
    res <- list(n0=out$n0,ntotal=out$ntotal,jj=out$jj,zl=out$zl,zr=out$zr,y=out$y,
                enc=out$enc,day=out$day,sigma=out$sigma,covar=out$covar,q=out$q,
                del=out$del,WAIC=out$WAIC,DIC=out$DIC,Dbar=out$Dbar,pd=out$pd,trace1=out$trace1,
                trace2=out$trace2,trace3=out$trace3,trace4=out$trace4,venc=out$venc,
                vday=out$vday,vsigma=out$vsigma,vcovar=out$vcovar,vdel=out$vdel,
                vq=out$vq,asr=out$asr,vasr=out$vasr,casr=out$casr,vcasr=out$vcasr,
                sr=out$sr1,vsr=out$sr,tha=out$tha,mutha=out$mutha,vartha=out$vartha,
                trace5=out$trace5)
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
                asr=out$asr,vasr=out$vasr,casr=out$casr,vcasr=out$vcasr,sr=out$sr,vsr=out$sr)
  }else{
    stop("Inputs are not consistent")
  }}
  return(res)
}

