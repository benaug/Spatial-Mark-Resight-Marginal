e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.OneStage <- function(data,inits=NA,M=NA){
  library(abind)
  #extract observed data
  y.mID <- data$y.mID #marked detections
  y.mnoID <- data$y.mnoID #marked with no ID samples
  y.all <- data$y.all #all samples
  n.marked <- data$n.marked
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- data$K
  K1D <- data$K1D
  locs <- data$locs
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  lam0 <- inits$lam0
  sigma <- inits$sigma
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  #but update s.inits for marked individuals before assigning latent detections
  idx <- which(rowSums(y.mID)>0)
  for(i in idx){
    trps <- matrix(X[which(y.mID[i,]>0),1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  #update using telemetry if you have it
  if(!is.null(dim(data$locs))){
    max.locs <- dim(locs)[2]
    if(n.marked>1){
      tel.inds <- which(rowSums(is.na(locs[,,1]))<max.locs)
      n.locs.ind <- rowSums(!is.na(locs[,,1]))
    }else{
      tel.inds <- which(sum(is.na(locs[,,1]))<max.locs)
      n.locs.ind <- sum(!is.na(locs[,,1]))
    }
    print("using telemetry to initialize telemetered s. Remove from data if not using in the model.")
    #update s starts for telemetry guys
    for(i in tel.inds){
      if(n.locs.ind[i]>1){
        s.init[i,] <- colMeans(locs[i,1:n.locs.ind[i],])
      }else{
        s.init[i,] <- locs[i,1,]
      }
      #make sure s is in state space
      if(s.init[i,1]<xlim[1]){
        s.init[i,1] <- xlim[1] + 0.01
      }
      if(s.init[i,1]>xlim[2]){
        s.init[i,1] <- xlim[2] - 0.01
      }
      if(s.init[i,2]<ylim[1]){
        s.init[i,2] <- ylim[1] + 0.01
      }
      if(s.init[i,2]>ylim[2]){
        s.init[i,2] <- ylim[2] - 0.01
      }
    }
    n.locs.ind <- n.locs.ind[tel.inds]
  }else{
    tel.inds <- NA
    n.locs.ind <- NA
  }
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  
  y.true <- y.mID #add marked but unidentified to initialize s
  for(j in 1:J){
    #add marked no ID
    prob <- lamd[1:n.marked,j]
    prob <- prob/sum(prob)
    y.true[1:n.marked,j] <- y.true[1:n.marked,j] + rmultinom(1,y.mnoID[j],prob=prob)
  }
  
  #update s for marked individuals
  idx <- which(rowSums(y.true)>0)
  for(i in idx){
    trps <- matrix(X[y.true[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  
  #allocate unmarked samples to initializes s.inits for unmarked
  y.true2 <- matrix(0,M,J)
  for(j in 1:J){
    #add marked no ID
    prob <- lamd[,j]
    prob <- prob/sum(prob)
    y.true2[,j] <- y.true2[,j] + rmultinom(1,y.all[j],prob=prob)
  }
  
  idx <- which(rowSums(y.true2)>0)
  rem.idx <- which(idx<=n.marked)
  if(length(rem.idx)>0){
    idx <- idx[-rem.idx]
  }
  for(i in idx){
    trps <- matrix(X[y.true2[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  z.init <- 1*(rowSums(y.true2)>0)
  z.init[1:n.marked] <- 1
  
  #If using a habitat mask, move any s's initialized in non-habitat above to closest habitat
  e2dist  <-  function (x, y){
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
  }
  getCell  <-  function(s,res,cells){
    cells[trunc(s[1]/res)+1,trunc(s[2]/res)+1]
  }
  alldists <- e2dist(s.init,data$dSS)
  alldists[,data$InSS==0] <- Inf
  for(i in 1:M){
    this.cell <- data$cells[trunc(s.init[i,1]/data$res)+1,trunc(s.init[i,2]/data$res)+1]
    if(data$InSS[this.cell]==0){
      cands <- alldists[i,]
      new.cell <- which(alldists[i,]==min(alldists[i,]))
      s.init[i,] <- data$dSS[new.cell,]
    }
  }
  
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  
  
  #check starting logProbs
  #marked with ID obs
  logProb <- array(0,dim=c(n.marked,J))
  for(i in 1:n.marked){
    for(j in 1:J){
      logProb[i,j] <- dpois(y.mID[i,j],lamd[i,j]*data$K1D[j],log=TRUE)
    }
  }
  #marked no ID obs
  logProb <- rep(0,J)
  if(n.marked>1){
    lamd.mnoID <- colSums(lamd[1:n.marked,])
  }else{
    lamd.mnoID <- lamd[n.marked,]
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marked no ID observations.")
  #um obs
  logProb <- rep(0,J)
  lamd.um <- colSums(lamd[(n.marked+1):M,])
  for(j in 1:J){
    logProb[j] <- dpois(y.mnoID[j],lamd.mnoID[j]*data$K1D[j])
  }
 
  return(list(s=s.init,z=z.init,K1D=K1D,
              y.mID=y.mID,y.mnoID=y.mnoID,y.all=y.all,
              xlim=xlim,ylim=ylim,locs=locs,tel.inds=tel.inds,n.locs.ind=n.locs.ind))

}