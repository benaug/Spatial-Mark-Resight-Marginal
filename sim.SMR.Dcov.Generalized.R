e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SMR.Dcov.Generalized <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           lam0=NA,theta.d=NA,p0=NA,sigma=NA,K.mark=NA,K.sight=NA,
           X.mark=NA,X.sight=NA,xlim=NA,ylim=NA,K1D.mark=NA,K1D.sight=NA,
           theta.marked=c(1,0,0),theta.unmarked=1,tlocs=0,obstype="poisson"){
    if(sum(theta.marked)!=1)stop("theta.marked must sum to 1.")
    if(theta.unmarked<0|theta.unmarked>1)stop("theta.unmarked must be between 0 and 1.")
    library(abind)
    #get expected N
    cellArea <- res^2
    lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
    lambda.N <- sum(lambda.cell)
    #simulate realized N
    N <- rpois(1,lambda.N)
    
    #recreate some Dcov things so we can pass fewer arguments into this function
    x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
    y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
    dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
    cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
    n.cells <- nrow(dSS)
    n.cells.x <- length(x.vals)
    n.cells.y <- length(y.vals)
    
    # simulate a population of activity centers
    pi.cell <- lambda.cell/sum(lambda.cell)
    s.cell <- sample(1:n.cells,N,prob=pi.cell,replace=TRUE)
    #distribute activity centers uniformly inside cells
    s <- matrix(NA,nrow=N,ncol=2)
    for(i in 1:N){
      tmp <- which(cells==s.cell[i],arr.ind=TRUE) #x and y number
      s[i,1] <- runif(1,x.vals[tmp[1]]-res/2,x.vals[tmp[1]+res/2])
      s[i,2] <- runif(1,y.vals[tmp[2]]-res/2,y.vals[tmp[2]+res/2])
    }
    
    D.mark <- e2dist(s,X.mark)
    J.mark <- nrow(X.mark)
    D.sight <- e2dist(s,X.sight)
    J.sight <- nrow(X.sight)
    
    #trap operation - marking process
    if(!any(is.na(K1D.mark))){
      if(any(K1D.mark>K.mark)){
        stop("Some entries in K1D.mark are greater than K.mark.")
      }
      if(is.null(dim(K1D.mark))){
        if(length(K1D.mark)!=J.mark){
          stop("K1D.mark vector must be of length J.mark.")
        }
      }
    }else{
      print("K1D.mark not provided, assuming trap operation is perfect.")
      K1D.mark <- rep(K.mark,J.mark)
    }
    
    #trap operation - sighting process
    if(!any(is.na(K1D.sight))){
      if(any(K1D.sight>K.sight)){
        stop("Some entries in K1D.sight are greater than K.sight.")
      }
      if(is.null(dim(K1D.sight))){
        if(length(K1D.sight)!=J.sight){
          stop("K1D.sight vector must be of length J.sight.")
        }
      }
    }else{
      print("K1D.sight not provided, assuming trap operation is perfect.")
      K1D.sight <- rep(K.sight,J.sight)
    }
    
    #Capture and mark individuals
    y.mark <- array(0,dim=c(N,J.mark,K.mark))
    pd <- p0*exp(-D.mark*D.mark/(2*sigma*sigma))
    for(i in 1:N){
      for(j in 1:J.mark){
        y.mark[i,j,1:K1D.mark[j]] <- rbinom(K1D.mark[j],size=1,prob=pd[i,j])
      }
    }
    
    #resight individuals
    y <- array(0,dim=c(N,J.sight,K.sight))
    if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J.sight){
          y[i,j,1:K1D.sight[j]] <- rpois(K1D.sight[j],lamd[i,j])
        }
      } 
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J.sight){
          for(k in 1:K.sight){
            y[i,j,1:K1D.sight[j]] <- rnbinom(K1D.sight[j],mu=lamd[i,j],size=theta.d)
          }
        }
      } 
    }else{
      stop("obstype not recognized")
    }
    
    if(sum(y)==0)stop("No individuals captured. Reconsider parameter settings.")
    
    #get marked individuals
    cap.idx <- which(rowSums(y.mark)>0)
    if(length(cap.idx)==0)stop("Simulated 0 marked individuals.")
    
    y.mark <- y.mark[cap.idx,,] #marking process history
    n.marked <- length(cap.idx) #number of marked individuals
    #rearrange sighting history to put marked individuals at the top for code below to work correctly
    umguys <- setdiff(1:N,cap.idx) 
    y <- abind(y[cap.idx,,],y[umguys,,],along=1)
    #rearrange s to simulate telemetry correctly
    s <- rbind(s[cap.idx,],s[umguys,])
    
    #proceed with regular SMR simulation code
    IDmarked <- 1:n.marked
    
    y.event <- array(0,dim=c(N,J.sight,K.sight,3))
    #loop over cells with positive counts
    idx <- which(y>0,arr.ind=TRUE)
    for(l in 1:nrow(idx)){
      if(idx[l,1]<=n.marked){
        y.event[idx[l,1],idx[l,2],idx[l,3],] <- rmultinom(1,y[idx[l,1],idx[l,2],idx[l,3]],theta.marked)
      }else{
        y.event[idx[l,1],idx[l,2],idx[l,3],] <- rmultinom(1,y[idx[l,1],idx[l,2],idx[l,3]],c(0,theta.unmarked,1-theta.unmarked))
      }
    }
    
    y.mID <- apply(y.event[1:n.marked,,,1],c(1,2),sum)
    y.mnoID <- apply(y.event[1:n.marked,,,2],2,sum)
    y.um <- apply(y.event[(n.marked+1):N,,,2],2,sum)
    y.unk <- apply(y.event[1:n.marked,,,3],2,sum) + apply(y.event[(n.marked+1):N,,,3],2,sum)
    
    if(!sum(y)==(sum(y.mID)+sum(y.mnoID)+sum(y.um)+sum(y.unk)))stop("data simulator bug")
    
    #Telemetry observations
    if(tlocs>0){
      locs <- array(NA,dim=c(n.marked,tlocs,2))
      for(i in 1:n.marked){
        for(j in 1:tlocs){
          locs[i,j,] <- c(rnorm(1,s[IDmarked[i],1],sigma),rnorm(1,s[IDmarked[i],2],sigma))
        }
      }
    }else{
      locs <- NA
    }
    
    if(n.marked>1){
      n.M <- sum(rowSums(y[1:n.marked,,])>0)
    }else{
      if(sum(y[1,,])>0){
        n.M <- 1
      }else{
        n.M <- 0
      }
    }
    if(n.marked<N){
      n.UM <- sum(rowSums(y[(n.marked+1):N,,])>0)
    }else{
      n.UM <- 0
    }
    
    s.marked <- s[IDmarked,]
    y.mark <- apply(y.mark,c(1,2),sum)
    
    out <- list(y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk, #observed data
                n.marked=n.marked,locs=locs,n.M=n.M,n.UM=n.UM,
                y=y,s=s,s.marked=s.marked,#true data
                X.mark=X.mark,K.mark=K.mark,K1D.mark=K1D.mark,
                X.sight=X.sight,K.sight=K.sight,K1D.sight=K1D.sight,
                xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
                n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
                D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N)
    return(out)
  }