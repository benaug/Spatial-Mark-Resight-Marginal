#more complicated than necessary when using marginal approach, but gets the job done.

e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SMR.Dcov.Interspersed <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,n.marked=NA,lam0=NA,
           theta.d=NA,sigma=0.50,K=10,X=X,xlim=NA,ylim=NA,
           theta.marked=c(1,0,0),theta.unmarked=1,K2D=NA,tlocs=0,obstype="poisson",
           p.half=NA){
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
    
    D <- e2dist(s,X)
    J <- nrow(X)
    
    #trap operation
    if(any(is.na(K2D))){
      print("K2D not provided, assuming trap operation is perfect.")
      K2D <- matrix(1,J,K)
    }
    
    # Capture and mark individuals
    y <- array(0,dim=c(N,J,K))
    if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          y[i,j,1:K] <- rpois(K,lamd[i,j]*K2D[j,])
        }
      }
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          y[i,j,1:K] <- rnbinom(K,mu=lamd[i,j]*K2D[j,],size=theta.d)
        }
      } 
    }else{
      stop("obstype not recognized")
    }
    
    if(sum(y)==0)stop("No individuals captured. Reconsider parameter settings.")

    #simulate marked status
    marked.status <- matrix(1,n.marked,K)
    k.half <- floor(K/2)
    half.ind <- rep(0,n.marked)
    for(i in 1:n.marked){
      half.ind[i] <- rbinom(1,1,p.half) #1 is a "half individual" only marked for 1st half of occasions
      if(half.ind[i]==1){
        marked.status[i,1:k.half] <- 0
      }
    }    
    
    #can use data simulator that does not consider marked status from here
    #by moving marked.status=0 samples to new, fake, unmarked individuals
    n.half.inds <- sum(half.ind)
    y <- abind(y,array(0,dim=c(n.half.inds,J,K)),along=1)
    idx <- N + 1
    if(n.half.inds>0){
      these.inds <- which(half.ind==1)
      for(i in these.inds){
        tmp <- y[i,,1:k.half] #this marked guy's history to move to a new, fake, unmarked individual
        y[i,,1:k.half] <- 0 #zero marked guy out
        y[idx,,1:k.half] <- tmp #add to fake unmarked guy
        idx <- idx + 1
      }
    }
    N.original <- N #store this
    N <- N + n.half.inds
    
    IDmarked <- 1:n.marked
    
    y.event <- array(0,dim=c(N,J,K,3))
    #loop over cells with positive counts
    idx <- which(y>0,arr.ind=TRUE)
    for(l in 1:nrow(idx)){
      if(idx[l,1]<=n.marked){
        y.event[idx[l,1],idx[l,2],idx[l,3],] <- rmultinom(1,y[idx[l,1],idx[l,2],idx[l,3]],theta.marked)
      }else{
        y.event[idx[l,1],idx[l,2],idx[l,3],] <- rmultinom(1,y[idx[l,1],idx[l,2],idx[l,3]],c(0,theta.unmarked,1-theta.unmarked))
      }
    }
    
    y.mID <- y.event[1:n.marked,,,1]
    y.mnoID <- apply(y.event[1:n.marked,,,2],c(2,3),sum)
    y.um <- apply(y.event[(n.marked+1):N,,,2],c(2,3),sum)
    y.unk <- apply(y.event[1:n.marked,,,3],c(2,3),sum) + apply(y.event[(n.marked+1):N,,,3],c(2,3),sum)
    
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
    
    out <- list(y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,marked.status=marked.status, #observed data
                n.marked=n.marked,locs=locs,n.M=n.M,n.UM=n.UM,
                y=y,s=s,#true data
                X=X,K=K,K2D=K2D,
                xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
                n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
                D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N.original,lambda.N=lambda.N)
    return(out)
  }