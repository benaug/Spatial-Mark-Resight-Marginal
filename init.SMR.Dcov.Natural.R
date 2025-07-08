e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.Dcov.Natural <- function(data,inits=NA,M1=NA,M2=NA){
  library(abind)
  #extract observed data
  y.mID <- data$y.mID #marked detections
  y.mnoID <- data$y.mnoID #marked with no ID samples
  y.um <- data$y.um #unmarked samples
  y.unk <- data$y.unk #unknown marked status samples
  n.marked <- data$n.marked
  X <- as.matrix(data$X)
  J <- nrow(X)
  K <- data$K
  K1D <- data$K1D
  
  xlim <- data$xlim
  ylim <- data$ylim
  M.both <- M1 + M2
  
  #augment y.mID
  y.mID.aug <- matrix(0,M1,J)
  y.mID.aug[1:n.marked,] <- y.mID
  y.mID <- y.mID.aug
  ##pull out initial values
  lam0 <- inits$lam0
  sigma <- inits$sigma
  #assign random locations to assign latent ID samples to individuals
  s.init <- cbind(runif(M.both,xlim[1],xlim[2]), runif(M.both,ylim[1],ylim[2]))
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
  D <- e2dist(s.init, X)
  lamd <- lam0*exp(-D*D/(2*sigma*sigma))
  y.full <- matrix(0,M.both,J)
  y.full[1:M1,] <- y.mID
  for(j in 1:J){
    #add marked no ID
    prob <- lamd[1:M1,j]
    prob <- prob/sum(prob)
    y.full[1:M1,j] <- y.full[1:M1,j] + rmultinom(1,y.mnoID[j],prob=prob)
    #add unmarked
    prob <- c(rep(0,M1),lamd[(M1+1):M.both])
    prob <- prob/sum(prob)
    y.full[,j] <- y.full[,j] + rmultinom(1,y.um[j],prob=prob)
    #add unk
    prob <- lamd[,j]
    prob <- prob/sum(prob)
    y.full[,j] <- y.full[,j] + rmultinom(1,y.unk[j],prob=prob)
    
  }
  z.init <- 1*(rowSums(y.full)>0)
  z.init[1:n.marked] <- 1
  
  #update s for individuals assigned samples
  y.full2D <- y.full
  idx <- which(rowSums(y.full2D)>0)
  for(i in idx){
    trps <- matrix(X[y.full2D[i,]>0,1:2],ncol=2,byrow=FALSE)
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
  for(i in 1:M.both){
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
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marked with ID observations.")
  #marked no ID obs
  logProb <- rep(0,J)
  lamd.mnoID <- colSums(lamd[1:n.marked,])
  for(j in 1:J){
    logProb[j] <- dpois(y.mnoID[j],lamd.mnoID[j]*data$K1D[j])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Marked no ID observations.")
  #um obs
  logProb <- rep(0,J)
  lamd.um <- colSums(lamd[(M1+1):M.both,])
  for(j in 1:J){
    logProb[j] <- dpois(y.um[j],lamd.um[j]*data$K1D[j])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Unmarked observations.")
  #unk obs
  logProb <- rep(0,J)
  lamd.unk <- colSums(lamd[1:M.both,])
  for(j in 1:J){
    logProb[j] <- dpois(y.unk[j],lamd.unk[j]*data$K1D[j])
  }
  if(!is.finite(sum(logProb)))stop("Starting observation model likelihood not finite. Unknown marked status observations.")
  
  return(list(s=s.init,z=z.init,K1D=K1D,
              y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,xlim=xlim,ylim=ylim))

}