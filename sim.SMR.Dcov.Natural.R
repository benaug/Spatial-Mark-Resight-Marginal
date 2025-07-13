e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SMR.Dcov.Natural <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,p.marked=NA,lam0=NA,
           theta.d=NA,sigma=0.50,K=10,X=X,xlim=NA,ylim=NA,
           theta.marked=c(1,0,0),theta.unmarked=1,K1D=NA,obstype="poisson"){
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
    if(!any(is.na(K1D))){
      if(any(K1D>K)){
        stop("Some entries in K1D are greater than K.")
      }
      if(is.null(dim(K1D))){
        if(length(K1D)!=J){
          stop("K1D vector must be of length J.")
        }
      }
    }else{
      print("K1D not provided, assuming trap operation is perfect.")
      K1D <- rep(K,J)
    }
    
    # Capture and mark individuals
    y <- array(0,dim=c(N,J,K))
    if(obstype=="poisson"){
      if(is.na(lam0))stop("must provide lam0 for poisson obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          y[i,j,1:K1D[j]] <- rpois(K1D[j],lamd[i,j])
        }
      } 
    }else if(obstype=="negbin"){
      if(is.na(lam0))stop("must provide lam0 for negbin obstype")
      if(is.na(theta.d))stop("Must provide theta.d for negbin obstype")
      lamd <- lam0*exp(-D*D/(2*sigma*sigma))
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,1:K1D[j]] <- rnbinom(K1D[j],mu=lamd[i,j],size=theta.d)
          }
        }
      } 
    }else{
      stop("obstype not recognized")
    }
    
    if(sum(y)==0)stop("No individuals captured. Reconsider parameter settings.")
    
    
    cap.idx <- which(rowSums(y)>0)
    y <- y[cap.idx,,]
    #number of identifiable individuals who were captured. If we don't simulate a marked with ID detection for an
    #identifiable individual, we need to remove these from n.marked and y.mID (at the end of script)
    n.cap <- length(cap.idx)
    n.marked <- rbinom(1,n.cap,p.marked)
    IDmarked <- 1:n.marked
    
    y.event <- array(0,dim=c(n.cap,J,K,3))
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
    y.um <- apply(y.event[(n.marked+1):n.cap,,,2],2,sum)
    y.unk <- apply(y.event[1:n.marked,,,3],2,sum) + apply(y.event[(n.marked+1):n.cap,,,3],2,sum)
    
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
      n.UM <- sum(rowSums(y[(n.marked+1):nrow(y),,])>0)
    }else{
      n.UM <- 0
    }
    
    
    #fix this
    
    #final adjustment. We simulated the number of marked, identifiable, individuals, n.marked. But we may not
    #obtain a marked with ID detection of them, so remove them from y.mID and n.marked
    rem.idx <- which(rowSums(y.mID)==0)
    if(length(rem.idx)>0){
      y.mID <- y.mID[-rem.idx,]
    }
    n.marked <- nrow(y.mID)
    
    out <- list(y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,n.marked=n.marked, #observed data
                y=y,s=s,n.UM=n.UM, #true data
                X=X,K=K,K1D=K1D,
                xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
                n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
                D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N)
    return(out)
  }