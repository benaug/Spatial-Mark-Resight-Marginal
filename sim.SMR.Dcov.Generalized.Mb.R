e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

sim.SMR.Dcov.Generalized.Mb <-
  function(D.beta0=NA,D.beta1=NA,D.cov=NA,InSS=NA,res=NA,
           lam0=NA,theta.d=NA,p0.p=NA,p0.c=NA,sigma=NA,K.mark=NA,K.sight=NA,
           X.mark=NA,X.sight=NA,xlim=NA,ylim=NA,K2D.mark=NA,K1D.sight=NA,
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
    #trap operation - sighting process
    if(any(is.na(K2D.mark))){
      print("K2D.mark not provided, assuming trap operation is perfect.")
      K2D.mark <- matrix(1,J.mark,K.mark)
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
    y.state <- array(0,dim=c(N,J.mark,K.mark))
    pd.p <- p0.p*exp(-D.mark*D.mark/(2*sigma*sigma))
    pd.c <- p0.c*exp(-D.mark*D.mark/(2*sigma*sigma))
    for(i in 1:N){
      for(j in 1:J.mark){
        for(k in 1:K.mark){
          y.mark[i,j,k] <- rbinom(1,K2D.mark[j,k],pd.p[i,j]*(y.state[i,j,k]==0) + 
                               pd.c[i,j]*(y.state[i,j,k]==1))
          if(y.mark[i,j,k]==1&k!=K.mark){
            y.state[i,j,(k+1):K.mark] <- 1
          }
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
    umguys <- setdiff(1:N,IDmarked) 
    
    #split sightings into marked and unmarked histories
    y.marked <- y[IDmarked,,]
    if(length(IDmarked)==1){ #if only one marked guy, make y.marked an array again
      y.mark <- array(y.mark,dim=c(1,J.mark,K.mark))
      y.marked <- array(y.marked,dim=c(1,J.sight,K.sight))
    }
    n.samples <- sum(y[umguys,,])
    y.unmarked <- array(0,dim=c(n.samples,J.sight,K.sight))
    IDum <- rep(NA,n.samples)
    if(n.samples>0){
      idx <- 1
      for(i in 1:length(umguys)){
        for(j in 1:J.sight){ #then traps
          for(k in 1:K.sight){ #then occasions
            if(y[umguys[i],j,k]>0){ #is there at least one sample here?
              for(l in 1:y[umguys[i],j,k]){ #then samples
                y.unmarked[idx,j,k] <- 1
                IDum[idx] <- umguys[i]
                idx <- idx+1
              }
            }
          }
        }
      }
    }
    
    #ID/marked status observation model for marked individuals
    # if(theta.marked[1]!=0&sum(y.marked)>0){#bug fix from Glenn Stauffer
    if(sum(y.marked)>0){#bug fix from Glenn Stauffer. BA change 4/22/23. 
      idx1 <- which(y.marked>0)#used to extract counts
      count <- y.marked[idx1]
      idx2 <- which(y.marked>0,arr.ind=TRUE)#used to move counts
      idx3 <- rep(1:nrow(idx2),count)#repeat array indices for counts>1
      idx2 <- idx2[idx3,]
      if(!is.matrix(idx2)){
        idx2 <- matrix(idx2,ncol=3)
      }
      n.sample.mark <- sum(count)
      outcome <- rmultinom(n.sample.mark,1,theta.marked)
      mnoID.idx <- which(outcome[2,]==1)
      munk.idx <- which(outcome[3,]==1)
      
      #remove marked but no ID and unk marked status from y.marked
      #create y.marked.noID and y.unk
      nmnoID <- length(mnoID.idx)
      nmunk <- length(munk.idx)
      if(nmnoID>0){
        y.marked.noID <- array(0,dim=c(nmnoID,J.sight,K.sight))
        for(i in 1:length(mnoID.idx)){
          #delete this sighting
          y.marked[idx2[mnoID.idx[i],1],idx2[mnoID.idx[i],2],idx2[mnoID.idx[i],3]] <- 
            y.marked[idx2[mnoID.idx[i],1],idx2[mnoID.idx[i],2],idx2[mnoID.idx[i],3]]-1
          #add this sighting
          y.marked.noID[i,idx2[mnoID.idx[i],2],idx2[mnoID.idx[i],3]] <- 1
        }
        IDmnoID <- idx2[mnoID.idx,1]
      }else{
        y.marked.noID <- IDnoID <- NA
        IDmnoID <- NA
      }
      if(nmunk>0){
        y.unk <- array(0,dim=c(nmunk,J.sight,K.sight))
        for(i in 1:length(munk.idx)){
          #delete this sighting
          y.marked[idx2[munk.idx[i],1],idx2[munk.idx[i],2],idx2[munk.idx[i],3]] <- 
            y.marked[idx2[munk.idx[i],1],idx2[munk.idx[i],2],idx2[munk.idx[i],3]]-1
          #add this sighting
          y.unk[i,idx2[munk.idx[i],2],idx2[munk.idx[i],3]] <- 1
        }
        IDunk <- idx2[munk.idx,1]
        IDunkType <- rep("marked",length(munk.idx))
      }else{
        y.unk <- IDunk <- IDunkType <- NA
      }
    }else{
      y.marked.noID <- IDmnoID <- NA
      y.unk <- IDunk <- IDunkType <- NA
    }
    
    #marked status observation model for unmarked individuals
    if(theta.unmarked!=1){
      outcome <- rbinom(n.samples,1,theta.unmarked)
      unk.idx <- which(outcome==0)
      nunk <- length(unk.idx)
      if(nunk>0){
        #extract um history to unk
        y.unk2 <- y.unmarked[unk.idx,,]
        if(nunk==1){ #reformat to array
          y.unk2 <- array(y.unk2,dim=c(1,J.sight,K.sight))
        }
        IDunk2 <- IDum[unk.idx]
        IDunkType2 <- rep("unmarked",length(IDunk2))
        #remove unk from um history
        y.unmarked <- y.unmarked[-unk.idx,,]
        if(length(dim(y.unmarked))==2){ #reformat to array
          y.unmarked <- array(y.unmarked,dim=c(1,J.sight,K.sight))
        }
        IDum <- IDum[-unk.idx]
      }else{
        IDunk2 <- IDunkType2=NA
      }
    }else{
      y.unk2 <- NA
      IDunk2 <- IDunkType2 <- NA
    }
    
    #combine y.unk if there are members from both marked and unmarked
    if(!any(is.na(IDunk))&!any(is.na(IDunk2))){
      IDunk <- c(IDunk,IDunk2)
      IDunkType <- c(IDunkType,IDunkType2)
      y.unk <- abind(y.unk,y.unk2,along=1)
    }else if(!any(is.na(IDunk2))){#or if just unk from unmarked, rename them
      IDunk <- IDunk2
      IDunkType <- IDunkType2
      y.unk <- y.unk2
    }
    
    #check data
    y.check <- y*0
    for(i in 1:length(IDmarked)){
      y.check[IDmarked[i],,] <- y.marked[i,,]
    }
    if(length(IDum)>0){
      for(i in 1:length(IDum)){
        y.check[IDum[i],,] <- y.check[IDum[i],,]+y.unmarked[i,,]
      }
    }
    if(all(!is.na(IDunk))){
      for(i in 1:length(IDunk)){
        y.check[IDunk[i],,] <- y.check[IDunk[i],,]+y.unk[i,,]
      }
    }
    if(theta.marked[2]>0){
      if(all(is.finite(IDmnoID))){
        for(i in 1:length(IDmnoID)){
          y.check[IDmnoID[i],,] <- y.check[IDmnoID[i],,]+y.marked.noID[i,,]
        }
      }
    }
    if(!all(y==y.check)){
      stop("Error rebuilding data. Report bug.")
    }
    dimnames(y.unk) <- NULL
    
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
    
    #bug check add locs
    # for(i in 1:n.marked){
    #   points(locs[i,,1],locs[i,,2],pch=".")
    # }
    
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
    
    #formatting for marginal sampler
    y.marked <- apply(y.marked,c(1,2),sum)
    y.um <- apply(y.unmarked,2,sum)
    if(!all(is.na(y.marked.noID))){
      y.mnoID <- apply(y.marked.noID,2,sum)
    }else{
      y.mnoID <- rep(0,J.sight)
    }
    if(!all(is.na(y.unk))){
      y.unk <- apply(y.unk,2,sum)
    }else{
      y.unk <- rep(0,J.sight)
    }
    s.marked <- s[IDmarked,]
    # y.mark <- apply(y.mark,c(1,2),sum) #leave 3D
    
    out <- list(y.mark=y.mark,y.mID=y.marked,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk, #observed data
                n.marked=n.marked,locs=locs,n.M=n.M,n.UM=n.UM,
                y=y,s=s,s.marked=s.marked,#true data
                X.mark=X.mark,K.mark=K.mark,K2D.mark=K2D.mark,
                X.sight=X.sight,K.sight=K.sight,K1D.sight=K1D.sight,
                xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
                n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
                D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N)
    return(out)
  }