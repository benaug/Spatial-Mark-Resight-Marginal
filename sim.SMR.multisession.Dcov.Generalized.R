e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SMR.multisession.Dcov.Generalized <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           lam0=NA,p0=NA,theta.d=NA,sigma=NA,K.mark=NA,K.sight=NA,
           X.mark=X.mark,X.sight=X.sight,xlim=NA,ylim=NA,
           theta.marked=NA,theta.unmarked=NA,K1D.mark=NA,K1D.sight=NA,tlocs=NA,obstype="poisson"){
    if(length(D.beta0)!=N.session)stop("D.beta0 must be of length N.session")
    if(length(D.beta1)!=N.session)stop("D.beta1 must be of length N.session")
    if(length(lam0)!=N.session)stop("lam0 must be of length N.session")
    if(length(p0)!=N.session)stop("n.marked must be of length N.session")
    if(length(sigma)!=N.session)stop("sigma must be of length N.session")
    if(length(K.mark)!=N.session)stop("K.mark must be of length N.session")
    if(length(K.sight)!=N.session)stop("K.sight must be of length N.session")
    if(length(X.mark)!=N.session)stop("X.mark must be of length N.session")
    if(length(X.sight)!=N.session)stop("X.sight must be of length N.session")
    if(!is.matrix(theta.marked))stop("theta.marked must be a matrix")
    if(nrow(theta.marked)!=N.session)stop("theta.marked must have N.session rows")
    if(ncol(theta.marked)!=3)stop("theta.marked must have N.session columns")
    if(!all(rowSums(theta.marked)==1))stop("theta.marked rows must all sum to 1.")
    if(length(theta.unmarked)!=N.session)stop("theta.unmarked must be of length N.session")
    if(obstype=="negbin"){
      if(!any(is.na(theta.d))){
        if(length(theta.d)!=N.session)stop("theta.d must be of length N.session")
      }else{
        theta.d <- rep(NA,N.session)
      }
    }else{
      theta.d <- rep(NA,N.session)
    }
    if(length(theta.unmarked)!=N.session)stop("theta.unmarked must be of length N.session")
    
    J.mark <- J.sight <- rep(NA,N.session)
    for(g in 1:N.session){
      X.mark[[g]] <- as.matrix(X.mark[[g]])
      X.sight[[g]] <- as.matrix(X.sight[[g]])
      J.mark[g] <- nrow(X.mark[[g]])
      J.sight[g] <- nrow(X.sight[[g]])
    }
    
    if(any(is.na(K1D.mark))){
      print("K1D.mark not provided, assuming trap operation is perfect.")
      K1D.mark <- vector("list",N.session)
      for(g in 1:N.session){
        K1D.mark[[g]] <- rep(K.mark[g],J.mark[g])
      }
    }
    if(any(is.na(K1D.sight))){
      print("K1D.sight not provided, assuming trap operation is perfect.")
      K1D.sight <- vector("list",N.session)
      for(g in 1:N.session){
        K1D.sight[[g]] <- rep(K.sight[g],J.sight[g])
      }
    }
    
    #simulate sessions one at a time
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <- sim.SMR.Dcov.Generalized(D.beta0=D.beta0[g],D.beta1=D.beta1[g],res=res[g],D.cov=D.cov[[g]],InSS=InSS[[g]],
                                            theta.marked=theta.marked[g,],theta.unmarked=theta.unmarked[g],
                                            lam0=lam0[g],p0=p0[g],sigma=sigma[g],theta.d=theta.d[g],K.mark=K.mark[g],K.sight=K.sight[g],
                                            K1D.mark=K1D.mark[[g]],K1D.sight=K1D.sight[[g]],
                                            X.mark=X.mark[[g]],X.sight=X.sight[[g]],xlim=xlim[g,],ylim=ylim[g,],tlocs=tlocs[g],
                                            obstype=obstype)
    }
    return(data)
  }