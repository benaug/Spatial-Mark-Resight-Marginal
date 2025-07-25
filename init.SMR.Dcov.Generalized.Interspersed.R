e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.Generalized.Interspersed <- function(data,inits=NA,M=NA){
  library(abind)
  #extract observed data
  y.mark <- data$y.mark
  y.mID <- data$y.mID #marked detections
  y.mnoID <- data$y.mnoID #marked with no ID samples
  y.um <- data$y.um #unmarked samples
  y.unk <- data$y.unk #unknown marked status samples
  marked.status <- data$marked.status
  n.marked <- data$n.marked
  X.mark <- as.matrix(data$X.mark)
  J.mark <- nrow(X.mark)
  K.mark <- data$K.mark
  K1D.mark <- data$K1D.mark
  X.sight <- as.matrix(data$X.sight)
  J.sight <- nrow(X.sight)
  K.sight <- data$K.sight
  K2D.sight <- data$K2D.sight
  locs <- data$locs
  
  xlim <- data$xlim
  ylim <- data$ylim
  
  ##pull out initial values
  p0 <- inits$p0
  lam0 <- inits$lam0
  sigma <- inits$sigma
  
  #augment marking process capture history
  y.aug <- matrix(0,M,J.mark)
  y.aug[1:n.marked,] <- y.mark
  y.mark <- y.aug
  
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M,xlim[1],xlim[2]), runif(M,ylim[1],ylim[2]))
  #but update s.inits for marked individuals before assigning latent detections
  y.mark2D <- apply(y.mark,c(1,2),sum)
  y.mID2D <-  apply(y.mID,c(1,2),sum)
  y.both <- cbind(y.mark2D[1:n.marked,],y.mID2D)
  X.both <- rbind(X.mark,X.sight)
  idx <- which(rowSums(y.both)>0)
  for(i in idx){
    trps <- matrix(X.both[which(y.both[i,]>0),1:2],ncol=2,byrow=FALSE)
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
    #update using telemetry if you have it
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
  
  #build plausible true sighting history to better initialize s
  D.sight <- e2dist(s.init, X.sight)
  lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
  y.true <- array(0,dim=c(M,J.sight,K.sight))
  y.true[1:n.marked,,] <- y.mID
  marked.status.full <- matrix(0,M,K.sight)
  marked.status.full[1:n.marked,] <- marked.status
  for(j in 1:J.sight){
    for(k in 1:K.sight){
      #add marked no ID
      prob <- lamd[1:n.marked,j]*marked.status[,k]
      prob <- prob/sum(prob)
      if(sum(y.true[1:n.marked,j,k])>0){
        y.true[1:n.marked,j,k] <- y.true[1:n.marked,j,k] + rmultinom(1,y.mnoID[j,k],prob=prob)
      }
      #add unmarked
      prob <- c(lamd[1:n.marked,j]*(1-marked.status[,k]),lamd[(n.marked+1):M,j])
      prob <- prob/sum(prob)
      y.true[,j,k] <- y.true[,j,k] + rmultinom(1,y.um[j,k],prob=prob)
      #add unk
      prob <- lamd[,j]
      prob <- prob/sum(prob)
      y.true[,j,k] <- y.true[,j,k] + rmultinom(1,y.unk[j,k],prob=prob)
    }
  }
  z.init <- 1*(rowSums(y.true)>0)
  z.init[1:n.marked] <- 1
  
  #update s.init given marking and sighting histories
  y.true2D <- apply(y.true,c(1,2),sum)
  y.both <- cbind(y.mark,y.true2D)
  X.both <- rbind(X.mark,X.sight)
  idx <- which(rowSums(y.both)>0)
  for(i in idx){
    trps <- matrix(X.both[y.both[i,]>0,1:2],ncol=2,byrow=FALSE)
    if(nrow(trps)>1){
      s.init[i,] <- c(mean(trps[,1]),mean(trps[,2]))
    }else{
      s.init[i,] <- trps
    }
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
  
  
  #check starting logProbs
  D.mark <- e2dist(s.init, X.mark)
  pd <- p0*exp(-D.mark*D.mark/(2*sigma*sigma))
  D.sight <- e2dist(s.init, X.sight)
  lamd <- lam0*exp(-D.sight*D.sight/(2*sigma*sigma))
  
  #marking process
  logProb <- array(0,dim=c(M,J.mark))
  for(i in 1:M){
    for(j in 1:J.mark){
      logProb[i,j] <- dbinom(y.mark[i,j],K1D.mark[j],pd[i,j],log=TRUE)
    }
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marking process.")
  
  #marked with ID obs
  logProb <- array(0,dim=c(n.marked,J.sight,K.sight))
  for(i in 1:n.marked){
    for(j in 1:J.sight){
      logProb[i,j,] <- dpois(y.mID[i,j,],lamd[i,j]*data$K2D.sight[j,]*inits$theta.marked[1],log=TRUE)
    }
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marked with ID observations.")
  #marked no ID obs
  logProb <- matrix(0,J.sight,K.sight)
  if(n.marked>1){
    lamd.mnoID <- colSums(lamd[1:n.marked,])
  }else{
    lamd.mnoID <- lamd[n.marked,]
  }
  for(j in 1:J.sight){
    logProb[j,] <- dpois(y.mnoID[j,],lamd.mnoID[j]*data$K2D.sight[j,]*inits$theta.marked[2])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marked no ID observations.")
  #um obs
  logProb <- matrix(0,J.sight,K.sight)
  lamd.um <- colSums(lamd[(n.marked+1):M,])
  for(j in 1:J.sight){
    logProb[j,] <- dpois(y.um[j,],lamd.um[j]*data$K2D.sight[j,]*inits$theta.unmarked[2])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Unmarked observations.")
  #unk obs
  logProb <- matrix(0,J.sight,K.sight)
  lamd.unk <- colSums(rbind(lamd[1:n.marked,]*inits$theta.marked[3],lamd[(n.marked+1):M,]*inits$theta.unmarked[3]))
  for(j in 1:J.sight){
    logProb[j,] <- dpois(y.unk[j,],lamd.unk[j]*data$K2D.sight[j,])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Unknown marked status observations.")
  
  return(list(s=s.init,z=z.init,K1D.mark=K1D.mark,K2D.sight=K2D.sight,
              y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,marked.status=data$marked.status,
              xlim=xlim,ylim=ylim,locs=locs,tel.inds=tel.inds,n.locs.ind=n.locs.ind))

}