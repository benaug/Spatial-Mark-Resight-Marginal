e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SMR.Dcov.Generalized.Interspersed.GPSfail <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           lam0=NA,theta.d=NA,p0=NA,sigma=NA,K.mark=NA,K.sight=NA,
           X.mark=NA,X.sight=NA,xlim=NA,ylim=NA,K1D.mark=NA,K2D.sight=NA,
           theta.marked1=c(1,0,0),theta.marked2=c(0,1,0),theta.unmarked=1,
           p.fail=NA,tlocs=0,obstype="poisson",
           K.order=NA){
    if(sum(theta.marked1)!=1)stop("theta.marked1 must sum to 1.")
    if(sum(theta.marked2)!=1)stop("theta.marked2 must sum to 1.")
    if(theta.marked2[1]!=0)stop("We are assuming you can't ID inds with dead GPS collars")
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
      s[i,1] <- runif(1,x.vals[tmp[1]]-res/2,x.vals[tmp[1]]+res/2)
      s[i,2] <- runif(1,y.vals[tmp[2]]-res/2,y.vals[tmp[2]]+res/2)
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
    if(any(is.na(K2D.sight))){
      print("K2D.sight not provided, assuming trap operation is perfect.")
      K2D.sight <- matrix(1,J.sight,K.sight)
    }
    
    #check capture order
    if(is.na(K.order[1])){
      warning("Assuming all marking occasions preceded the first sighting occasion since 'K.order' was omitted")
      K.order=c(rep("M",K.mark),rep("S",K.sight))
    }
    if(!all(c("M","S")%in%names(table(K.order)))|!all(names(table(K.order))%in%c("M","S"))){
      stop("K.order must only contain characters'M'and 'S' indicating marking and sighting sessions")
    }
    if(length(K.order)!=sum(K.mark,K.sight)){
      stop("K.order is not the right size")
    }
    if((sum(K.order=="M")!=K.mark)|(sum(K.order=="S")!=K.sight)){
      stop("Fix number of M and S in K.order to match K.mark and K.sight")
    }
    
    #Capture and mark individuals
    y.mark <- array(0,dim=c(N,J.mark,K.mark))
    pd <- p0*exp(-D.mark*D.mark/(2*sigma*sigma))
    for(i in 1:N){
      for(j in 1:J.mark){
        if(K1D.mark[j]>0){
          y.mark[i,j,1:K1D.mark[j]] <- rbinom(K1D.mark[j],size=1,prob=pd[i,j])
        }
      }
    }
    
    #resight individuals
    y <- array(0,dim=c(N,J.sight,K.sight))
    if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J.sight){
          y[i,j,1:K.sight] <- rpois(K.sight,lamd[i,j]*K2D.sight[j,])
        }
      } 
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J.sight){
          y[i,j,1:K.sight] <- rnbinom(K.sight,mu=lamd[i,j]*K2D.sight[j,],size=theta.d)
        }
      } 
    }else{
      stop("obstype not recognized")
    }
    
    if(sum(y)==0)stop("No individuals captured. Reconsider parameter settings.")
    
    #get marked individuals
    cap.idx <- which(rowSums(y.mark)>0)
    if(length(cap.idx)==0)stop("Simulated 0 marked individuals.")
    y.mark <- y.mark[cap.idx,,,drop=FALSE] #marking process history
    n.marked <- length(cap.idx) #number of marked individuals

    #rearrange sighting history to put marked individuals at the top for code below to work correctly
    umguys <- setdiff(1:N,cap.idx) 
    y <- abind(y[cap.idx,,,drop=FALSE],y[umguys,,,drop=FALSE],along=1)
    #rearrange s to simulate telemetry correctly
    s <- rbind(s[cap.idx,],s[umguys,])
    
    #Figure out GPS status on each sighting occasion
    #marked.status: 0 = not yet marked, 1 = working GPS, 2 = dead GPS
    markedT <- 1*(apply(y.mark,c(1,3),sum)>0)
    marked.status <- matrix(0,ncol=K.sight,nrow=n.marked)
    markocc <- which(K.order=="M")
    sightocc <- which(K.order=="S")
    for(i in 1:n.marked){
      gps.status <- 0
      for(k in 1:K.sight){
        if(k==1){
          newmark <- which(markocc<sightocc[k])
        }else{
          newmark <- which(markocc>sightocc[k-1] & markocc<sightocc[k])
        }
        #if this individual was captured/re-collared since the previous sighting occasion,
        #GPS is reset to working
        if(length(newmark)>0 && any(markedT[i,newmark]==1)){
          gps.status <- 1
          new.mark <- TRUE
        }else{
          new.mark <- FALSE
        }
        #if already marked and not newly marked, GPS can fail
        if(gps.status==1 && !new.mark){
          gps.status <- sample(c(1,2),1,prob=c(1-p.fail,p.fail))
        }
        marked.status[i,k] <- gps.status
      }
    }
   
    #can use data simulator that does not consider marked status from here
    #by moving marked.status=0 samples to new, fake, unmarked individuals
    partial.inds <- which(rowSums(marked.status!=0)<K.sight)
    n.partial.inds <- length(partial.inds)
    y.orig <- y
    y <- abind(y,array(0,dim=c(n.partial.inds,J.sight,K.sight)),along=1)
    idx <- N + 1
    if(n.partial.inds>0){
      for(i in partial.inds){
        k.split <- which(marked.status[i,]!=0)[1] - 1
        if(k.split>0){
          tmp <- y[i,,1:k.split] #this marked guy's pre-marking history to move to a new, fake, unmarked individual
          y[i,,1:k.split] <- 0 #zero pre-marking history out for marked guy
          y[idx,,1:k.split] <- tmp #add pre-marking history to fake unmarked guy
        }
        idx <- idx + 1
      }
    }
    
    N.original <- N #store this
    N <- N + n.partial.inds
    
    
    IDmarked <- 1:n.marked
    
    y.event <- array(0,dim=c(N,J.sight,K.sight,3))
    #loop over cells with positive counts
    idx <- which(y>0,arr.ind=TRUE)
    for(l in 1:nrow(idx)){
      i <- idx[l,1]
      j <- idx[l,2]
      k <- idx[l,3]
      if(i<=n.marked){
        if(marked.status[i,k]==1){
          y.event[i,j,k,] <- rmultinom(1,y[i,j,k],theta.marked1)
        }else if(marked.status[i,k]==2){
          y.event[i,j,k,] <- rmultinom(1,y[i,j,k],theta.marked2)
        }else{
          stop("Simulator bug: marked individual has positive y while marked.status is 0.")
        }
      }else{
        y.event[i,j,k,] <- rmultinom(1,y[i,j,k],c(0,theta.unmarked,1-theta.unmarked))
      }
    }
    
    y.mID <- array(y.event[1:n.marked,,,1],dim=c(n.marked,J.sight,K.sight))
    y.mnoID <- apply(y.event[1:n.marked,,,2,drop=FALSE],c(2,3),sum)
    y.unk.marked <- apply(y.event[1:n.marked,,,3,drop=FALSE],c(2,3),sum)
    if(n.marked<N){
      y.um <- apply(y.event[(n.marked+1):N,,,2,drop=FALSE],c(2,3),sum)
      y.unk.unmarked <- apply(y.event[(n.marked+1):N,,,3,drop=FALSE],c(2,3),sum)
    }else{
      y.um <- matrix(0,J.sight,K.sight)
      y.unk.unmarked <- matrix(0,J.sight,K.sight)
    }
    y.unk <- y.unk.marked + y.unk.unmarked
    
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
    
    n.M <- sum(apply(y.orig[1:n.marked,,,drop=FALSE],1,sum)>0)
    if(n.marked<N.original){
      n.UM <- sum(apply(y.orig[(n.marked+1):N.original,,,drop=FALSE],1,sum)>0)
    }else{
      n.UM <- 0
    }
    
    y.mark <- apply(y.mark,c(1,2),sum)
    s.marked <- s[IDmarked,]
    
    out <- list(y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,marked.status=marked.status, #observed data
                n.marked=n.marked,locs=locs,n.M=n.M,n.UM=n.UM,
                y=y.orig,s=s,s.marked=s.marked,#true data
                X.mark=X.mark,K.mark=K.mark,K1D.mark=K1D.mark,
                X.sight=X.sight,K.sight=K.sight,K2D.sight=K2D.sight,
                xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
                n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
                D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N.original,lambda.N=lambda.N)
    return(out)
  }