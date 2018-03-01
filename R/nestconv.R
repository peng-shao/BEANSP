#' @title Data Format Conversion
#' @description Convert date format data or interval format data to the standard format.
#' @param type If data are interval format, type=0; if data are date format, type=1.
#' @param nn The number of observed nests (sample size).
#' @param jj The number of time units(days) a nest is required to survive to be considered successful.
#' @param data Input data.
#' @param censored If data are censored, censored=TRUE; if data are not censored, censored=FALSE
#' @return The returned values are \code{zl,zr}. If censored=FALSE, \code{ul,ur} are also returned.
#' @author
#' Chong He, Yiqun Yang, Jing Cao, Peng Shao
#' @keywords methods
#' @references Cao, J., He, C., Suedkamp Wells, K.M., Millspaugh, J.J., and
#' Ryan, M.R. (2009).  Modeling age and nest-specific survival using a hierarchical
#' Bayesian approach. Biometrics, 65, 1052-1062.
#' @export


nestconv <-
  function(type,nn,jj,data, censored=TRUE){
    if (censored){
      y<-rep(0,nn)
      zl<-rep(0,nn)
      zr<-rep(0,nn)
      if(type==0){
        # interval format
        for(i in 1:nn){
          if(data[i,4]=="S")
          {y[i]=1}
          else if(data[i,4]=="F")
          {y[i]=0}
          else
          {y[i]=2}
          zl[i]=data[i,2]
          zr[i]=data[i,3]
        }
      }
      else{
        # date format
        for(i in 1:nn){
          if(data[i,5]=="S")
          {y[i]=1}
          else if(data[i,5]=="F")
          {y[i]=0}
          else
          {y[i]=2}
          ### adjust cases with date2==date3
          if(as.numeric(as.Date(data[i,3],"%m/%d/%Y"))==as.numeric(as.Date(data[i,4],"%m/%d/%Y")))
          {data[i,3]=data[i,2]}
          ### zl
          if(as.numeric(as.Date(data[i,2],"%m/%d/%Y"))==as.numeric(as.Date(data[i,4],"%m/%d/%Y")))
          {zl[i]=1}
          else if(as.numeric(as.Date(data[i,2],"%m/%d/%Y"))==as.numeric(as.Date(data[i,3],"%m/%d/%Y")))
          {zl[i]=1}
          else
          {zl[i]=2+as.numeric(as.Date(data[i,3],"%m/%d/%Y"))-as.numeric(as.Date(data[i,2],"%m/%d/%Y"))}
          ### zr
          zr[i]=1+as.numeric(as.Date(data[i,4],"%m/%d/%Y"))-as.numeric(as.Date(data[i,2],"%m/%d/%Y"))
          #### adjust for different jj
          if(zr[i]>jj) {
            zr[i]=jj
            if(zl[i]>jj) {zl[i]=jj}
          }
        }
      }
      id<-data[,1]
      temp<-cbind.data.frame(id,zl,zr,y)
    }else{
      y<-rep(0,nn)
      ul<-rep(0,nn)
      ur<-rep(0,nn)
      zl<-rep(0,nn)
      zr<-rep(0,nn)
      if(type==0){
        # interval format
        for(i in 1:nn){
          if(data[i,2]=="F"){y[i]=0}
          else {y[i]=1}
          ul[i]=data[i,3]
          ur[i]=data[i,4]
          zl[i]=data[i,5]
          zr[i]=data[i,6]
        }
      }
      else{
        # date format
        for(i in 1:nn){
          if(data[i,5]=="F") {y[i]=0}
          else {y[i]=1}
          ### zl
          if(as.numeric(as.Date(data[i,2],"%m/%d/%Y"))==as.numeric(as.Date(data[i,4],"%m/%d/%Y")))
          {zl[i]=1}
          else
          {zl[i]=2+as.numeric(as.Date(data[i,3],"%m/%d/%Y"))-as.numeric(as.Date(data[i,2],"%m/%d/%Y"))}
          ### zr
          zr[i]=1+as.numeric(as.Date(data[i,4],"%m/%d/%Y"))-as.numeric(as.Date(data[i,2],"%m/%d/%Y"))
          #### adjust for different jj
          if(zr[i]>jj) {
            zr[i]=jj
            if(zl[i]>jj) {zl[i]=jj}
          }
          ### ul
          if(y[i]==0) {ul[i]=1}
          else {ul[i]=jj-zr[i]+1}
          ### ur
          ur[i]=jj-zl[i]+1
        }
      }
      id<-seq(1,nn)
      for(i in 1:nn){
        ### update the bounds
        if(ur[i]>(jj-zl[i]+1)){
          ur[i]=(jj-zl[i]+1)
          if(ul[i]>(jj-zl[i]+1)){
            ul[i]=jj-zl[i]+1
          }
        }
        if(y[i]==1){
          if(ul[i]<(jj-zr[i]+1)){
            ul[i]=jj-zr[i]+1
          }
        }
      }
      temp<-cbind.data.frame(id,ul,ur,zl,zr,y)
    }
    return(temp)
  }
