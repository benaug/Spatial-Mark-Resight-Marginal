e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}


sim.SMR.multisession.Dcov.Interspersed <-
  function(N.session=NA,D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           n.marked=NA,lam0=NA,theta.d=NA,sigma=NA,K=NA,X=X,xlim=NA,ylim=NA,
           theta.marked=NA,theta.unmarked=NA,K2D=NA,tlocs=NA,obstype="poisson",p.half=NA){
    if(length(D.beta0)!=N.session)stop("D.beta0 must be of length N.session")
    if(length(D.beta1)!=N.session)stop("D.beta1 must be of length N.session")
    if(length(n.marked)!=N.session)stop("n.marked must be of length N.session")
    if(length(lam0)!=N.session)stop("lam0 must be of length N.session")
    # if(length(theta.d)!=N.session)stop("theta.d must be of length N.session")
    if(length(sigma)!=N.session)stop("sigma must be of length N.session")
    if(length(K)!=N.session)stop("K must be of length N.session")
    if(length(X)!=N.session)stop("X must be of length N.session")
    if(!is.matrix(theta.marked))stop("theta.marked must be a matrix")
    if(nrow(theta.marked)!=N.session)stop("theta.marked must have N.session rows")
    if(ncol(theta.marked)!=3)stop("theta.marked must have N.session columns")
    if(!all(rowSums(theta.marked)==1))stop("theta.marked rows must all sum to 1.")
    if(length(theta.unmarked)!=N.session)stop("theta.unmarked must be of length N.session")
    if(obstype=="negbin"){
      if(length(theta.d)!=N.session)stop("theta.d must be of length N.session")
    }else{
      theta.d=rep(theta.d,N.session)
    }
    if(length(theta.unmarked)!=N.session)stop("theta.unmarked must be of length N.session")
    
    J <- rep(NA,N.session)
    for(g in 1:N.session){
      X[[g]] <- as.matrix(X[[g]])
      J[g] <- nrow(X[[g]])
    }
    
    #simulate sessions one at a time
    data <- vector("list",N.session)
    for(g in 1:N.session){
      data[[g]] <- sim.SMR.Dcov.Interspersed(D.beta0=D.beta0[g],D.beta1=D.beta1[g],D.cov=D.cov[[g]],InSS=InSS[[g]],res=res[g],
                                n.marked=n.marked[g],theta.marked=theta.marked[g,],theta.unmarked=theta.unmarked[g],
                                lam0=lam0[g],sigma=sigma[g],K=K[g],X=X[[g]],xlim=xlim[g,],ylim=ylim[g,],tlocs=tlocs[g],
                                obstype=obstype,theta.d=theta.d[g],p.half=p.half)
    }
    return(data)
  }