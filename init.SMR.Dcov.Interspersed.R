e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.Interspersed <- function(data,inits=NA,M=NA){
  library(abind)
  #extract observed data
  y.mID <- data$y.mID #marked detections
  y.mnoID <- data$y.mnoID #marked with no ID samples
  y.um <- data$y.um #unmarked samples
  y.unk <- data$y.unk #unknown marked status samples
  marked.status <- data$marked.status
  n.marked <- data$n.marked
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- data$K
  K2D <- data$K2D
  locs <- data$locs
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  lam0 <- inits$lam0
  sigma <- inits$sigma
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  y.full <- array(0,dim=c(M,J,K))
  y.full[1:n.marked,,] <- y.mID
  marked.status.full <- matrix(0,M,K)
  marked.status.full[1:n.marked,] <- marked.status
  for(j in 1:J){
    for(k in 1:K){
      #add marked no ID
      prob <- lamd[1:n.marked,j]*marked.status[,k]
      prob <- prob/sum(prob)
      if(sum(y.full[1:n.marked,j,k]>0)){
        y.full[1:n.marked,j,k] <- y.full[1:n.marked,j,k] + rmultinom(1,y.mnoID[j,k],prob=prob)
      }
      #add unmarked
      prob <- c(lamd[1:n.marked,j]*(1-marked.status[,k]),lamd[(n.marked+1):M])
      prob <- prob/sum(prob)
      y.full[,j,k] <- y.full[,j,k] + rmultinom(1,y.um[j,k],prob=prob)
      #add unk
      prob <- lamd[,j]
      prob <- prob/sum(prob)
      y.full[,j,k] <- y.full[,j,k] + rmultinom(1,y.unk[j,k],prob=prob)
    }
  }
  z.init <- 1*(rowSums(y.full)>0)
  z.init[1:n.marked] <- 1
  
  #update s for individuals assigned samples
  y.full2D <- apply(y.full,c(1,2),sum)
  idx <- which(rowSums(y.full2D)>0)
  for(i in idx){
    trps <- matrix(X[y.full2D[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
  }
  
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
  logProb <- array(0,dim=c(n.marked,J,K))
  for(i in 1:n.marked){
    for(j in 1:J){
      logProb[i,j,] <- dpois(y.mID[i,j,],lamd[i,j]*data$K2D[j,]*inits$theta.marked[1],log=TRUE)
    }
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marked with ID observations.")
  #marked no ID obs
  logProb <- matrix(0,J,K)
  if(n.marked>1){
    lamd.mnoID <- colSums(lamd[1:n.marked,])
  }else{
    lamd.mnoID <- lamd[n.marked,]
  }
  for(j in 1:J){
    logProb[j,] <- dpois(y.mnoID[j,],lamd.mnoID[j]*data$K2D[j,]*inits$theta.marked[2])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marked no ID observations.")
  #um obs
  logProb <- matrix(0,J,K)
  lamd.um <- colSums(lamd[(n.marked+1):M,])
  for(j in 1:J){
    logProb[j,] <- dpois(y.um[j,],lamd.um[j]*data$K2D[j,]*inits$theta.unmarked[2])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Unmarked observations.")
  #unk obs
  logProb <- matrix(0,J,K)
  lamd.unk <- colSums(rbind(lamd[1:n.marked,]*inits$theta.marked[3],lamd[(n.marked+1):M,]*inits$theta.unmarked[3]))
  for(j in 1:J){
    logProb[j,] <- dpois(y.unk[j,],lamd.unk[j]*data$K2D[j,])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Unknown marked status observations.")
  
  return(list(s=s.init,z=z.init,K2D=K2D,
              y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,marked.status=data$marked.status,
              xlim=xlim,ylim=ylim,locs=locs,tel.inds=tel.inds,n.locs.ind=n.locs.ind))
  
}