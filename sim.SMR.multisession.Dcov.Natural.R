e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


sim.SMR.multisession.Dcov.Natural <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           p.marked=NA,lam0=NA,theta.d=NA,sigma=NA,K=NA,X=X,xlim=NA,ylim=NA,
           theta.marked=NA,theta.unmarked=NA,K1D=NA,obstype="poisson"){
    if(length(D.beta0)!=N.session)stop("D.beta0 must be of length N.session")
    if(length(D.beta1)!=N.session)stop("D.beta1 must be of length N.session")
    if(length(p.marked)!=N.session)stop("p.marked must be of length N.session")
    if(length(lam0)!=N.session)stop("lam0 must be of length N.session")
    if(length(sigma)!=N.session)stop("sigma must be of length N.session")
    if(length(K)!=N.session)stop("K must be of length N.session")
    if(length(X)!=N.session)stop("X must be of length N.session")
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
    
    J <- rep(NA,N.session)
    for(g in 1:N.session){
      X[[g]] <- as.matrix(X[[g]])
      J[g] <- nrow(X[[g]])
    }
    
    if(any(is.na(K1D))){
      print("K1D not provided, assuming trap operation is perfect.")
      K1D <- vector("list",N.session)
      for(g in 1:N.session){
        K1D[[g]] <- rep(K[g],J[g])
      }
    }
    
    #simulate sessions one at a time
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <- sim.SMR.Dcov.Natural(D.beta0=D.beta0[g],D.beta1=D.beta1[g],D.cov=D.cov[[g]],InSS=InSS[[g]],res=res[g],
                                p.marked=p.marked[g],theta.marked=theta.marked[g,],theta.unmarked=theta.unmarked[g],
                                lam0=lam0[g],sigma=sigma[g],K=K[g],K1D=K1D[[g]],X=X[[g]],xlim=xlim[g,],ylim=ylim[g,],
                                obstype=obstype,theta.d=theta.d[g])
    }
    return(data)
  }