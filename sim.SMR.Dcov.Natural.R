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
    #identifiable individual, we need to remove these from n.marked and y.marked (at the end of script)
    n.marked <- rbinom(1,length(cap.idx),p.marked)
    IDmarked <- 1:n.marked
    umguys <- setdiff(1:nrow(y),IDmarked)    
    
    
    # cap.idx <- which(rowSums(y)>0)
    # n.marked <- rbinom(1,length(cap.idx),p.marked)
    # IDmarked <- cap.idx[1:n.marked]
    # umguys <- setdiff(1:N,IDmarked)
    
    #split sightings into marked and unmarked histories
    y.marked <- y[IDmarked,,]
    if(length(IDmarked)==1){ #if only one marked guy, make y.marked an array again
      y.marked <- array(y.marked,dim=c(1,J,K))
    }
    n.samples <- sum(y[umguys,,])
    y.unmarked <- array(0,dim=c(n.samples,J,K))
    IDum <- rep(NA,n.samples)
    if(n.samples>0){
      idx <- 1
      for(i in 1:length(umguys)){
        for(j in 1:J){ #then traps
          for(k in 1:K){ #then occasions
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
        y.marked.noID <- array(0,dim=c(nmnoID,J,K))
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
        y.unk <- array(0,dim=c(nmunk,J,K))
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
          y.unk2 <- array(y.unk2,dim=c(1,J,K))
        }
        IDunk2 <- IDum[unk.idx]
        IDunkType2 <- rep("unmarked",length(IDunk2))
        #remove unk from um history
        y.unmarked <- y.unmarked[-unk.idx,,]
        if(length(dim(y.unmarked))==2){ #reformat to array
          y.unmarked <- array(y.unmarked,dim=c(1,J,K))
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
    
    #formatting for marginal sampler
    y.marked <- apply(y.marked,c(1,2),sum)
    y.um <- apply(y.unmarked,2,sum)
    if(!all(is.na(y.marked.noID))){
      y.mnoID <- apply(y.marked.noID,2,sum)
    }else{
      y.mnoID <- rep(0,J)
    }
    if(!all(is.na(y.unk))){
      y.unk <- apply(y.unk,2,sum)
    }else{
      y.unk <- rep(0,J)
    }
    
    #final adjustment. We simulated the number of marked, identifiable, individuals, n.marked. But we may not
    #obtain a marked with ID detection of them, so remove them from y.marked and n.marked
    rem.idx <- which(rowSums(y.marked)==0)
    if(length(rem.idx)>0){
      y.marked <- y.marked[-rem.idx,]
    }
    n.marked <- nrow(y.marked)
    
    out <- list(y.mID=y.marked,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,n.marked=n.marked, #observed data
                y=y,s=s,n.UM=n.UM, #true data
                X=X,K=K,K1D=K1D,
                xlim=xlim,ylim=ylim,x.vals=x.vals,y.vals=y.vals,dSS=dSS,cells=cells,
                n.cells=n.cells,n.cells.x=n.cells.x,n.cells.y=n.cells.y,s.cell=s.cell,
                D.cov=D.cov,InSS=InSS,res=res,cellArea=cellArea,N=N,lambda.N=lambda.N)
    return(out)
  }