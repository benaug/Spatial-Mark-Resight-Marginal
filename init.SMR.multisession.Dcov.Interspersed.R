e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.multisession.Dcov.Interspersed <- function(data,inits=NA,M=NA){
  N.session <- length(data)
  n.marked <- sapply(data,function(x){x$n.marked})
  if(length(M)!=N.session)stop("Must supply an M for each session.")
  
  init.session <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #lam0 and sigma inits vary by session
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.session[[g]] <- init.SMR.Dcov.Interspersed(data[[g]],inits.use,M=M[g])
  }
  anyTelemetry <- FALSE
  locs.use <- vector("list",N.session)
  tlocs.sess.max <- rep(NA,N.session)
  for(g in 1:N.session){
    if(all(is.na(init.session[[g]]$locs))){
      locs.use <- NA
    }else{
      if(n.marked[g]>1){
        tlocs.sess.max[g] <- max(rowSums(!is.na(init.session[[g]]$locs[1:n.marked[g],,1])))
        locs.use[[g]] <- init.session[[g]]$locs[1:n.marked[g],1:tlocs.sess.max[g],1:2]
      }else{
        tlocs.sess.max[g] <- sum(!is.na(init.session[[g]]$locs[1:n.marked[g],,1]))
        locs.use[[g]] <- array(init.session[[g]]$locs[1,1:tlocs.sess.max[g],1:2],dim=dim(init.session[[g]]$locs))
      }
      anyTelemetry <- TRUE
    }
  }
  
  if(anyTelemetry){
    n.tel.inds <- unlist(lapply(locs.use,nrow))
    n.locs.ind <- matrix(NA,N.session,max(n.tel.inds))
    tel.inds <- matrix(NA,N.session,max(n.tel.inds))
  }else{
    tel.inds <- n.locs.ind <- n.tel.inds <- NA
  }
 
  J <- unlist(lapply(data,function(x){nrow(x$X)}))
  K <- unlist(lapply(data,function(x){x$K}))
  n.marked <- unlist(lapply(data,function(x){x$n.marked}))
  maxM <- max(M)
  s <- array(NA,dim=c(N.session,maxM,2))
  z <- matrix(NA,N.session,maxM)
  K2D <- array(NA,dim=c(N.session,max(J),max(K)))
  y.mID <- array(NA,dim=c(N.session,maxM,max(J),max(K)))
  y.mnoID <- array(NA,dim=c(N.session,max(J),max(K)))
  y.um <- array(NA,dim=c(N.session,max(J),max(K)))
  y.unk <- array(NA,dim=c(N.session,max(J),max(K)))
  marked.status <- array(NA,dim=c(N.session,max(n.marked),max(K)))
  
  n.cells <- unlist(lapply(data,function(x){x$n.cells}))
  n.cells.x <- unlist(lapply(data,function(x){x$n.cells.x}))
  n.cells.y <- unlist(lapply(data,function(x){x$n.cells.y}))
  n.cells.max <- max(n.cells)
  n.cells.x.max <- max(n.cells.x)
  n.cells.y.max <- max(n.cells.y)
  res <- unlist(lapply(data,function(x){x$res}))
  cellArea <- res^2
  xlim <- ylim <- matrix(NA,N.session,2)
  x.vals <- matrix(NA,N.session,n.cells.x.max)
  y.vals <- matrix(NA,N.session,n.cells.y.max)
  dSS <- array(NA,dim=c(N.session,n.cells.max,2))
  InSS <- array(0,dim=c(N.session,n.cells.max))
  D.cov <- array(NA,dim=c(N.session,n.cells.max))
  cells <- array(0,dim=c(N.session,n.cells.x.max,n.cells.y.max))
  
  for(g in 1:N.session){
    s[g,1:M[g],] <- init.session[[g]]$s
    z[g,1:M[g]] <- init.session[[g]]$z
    K2D[g,1:J[g],1:K[g]] <- init.session[[g]]$K2D
    y.mID[g,1:n.marked[g],1:J[g],1:K[g]] <- init.session[[g]]$y.mID
    y.mnoID[g,1:J[g],1:K[g]] <- init.session[[g]]$y.mnoID
    y.um[g,1:J[g],1:K[g]] <- init.session[[g]]$y.um
    y.unk[g,1:J[g],1:K[g]] <- init.session[[g]]$y.unk
    marked.status[g,1:n.marked[g],1:K[g]] <- init.session[[g]]$marked.status
    xlim[g,] <- data[[g]]$xlim
    ylim[g,] <- data[[g]]$ylim
    x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
    y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
    dSS[g,1:n.cells[g],] <- data[[g]]$dSS
    InSS[g,1:n.cells[g]] <- data[[g]]$InSS
    D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
    cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
    if(anyTelemetry){
      tel.inds[g,1:n.tel.inds[g]] <- init.session[[g]]$tel.inds
      n.locs.ind[g,1:n.tel.inds[g]] <- init.session[[g]]$n.locs.ind
    }
  }
  if(anyTelemetry){
    #reformat locs.use actually
    locs.use2 <- array(NA,dim=c(N.session,max(n.tel.inds,na.rm=TRUE),max(n.locs.ind,na.rm=TRUE),2))
    for(g in 1:N.session){
      locs.use2[g,1:nrow(locs.use[[g]]),1:ncol(locs.use[[g]]),] <- locs.use[[g]]
    }
    #remove unused telemetry dimensions if not all marked individuals telemetered 
    rem.idx <- which(colSums(is.na(n.locs.ind))==N.session)
    if(length(rem.idx)>0){
      n.locs.ind <- n.locs.ind[,-rem.idx]
      tel.inds <- tel.inds[,-rem.idx]
    }
  }else{
    locs.use2 <- NA
  }
  
  #put X in ragged array
  X.new <- array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],] <- data[[g]]$X
  }
  dummy.data <- matrix(0,N.session,maxM) #dummy data not used, doesn't really matter what the values are
  

  return(list(y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,marked.status=marked.status,
              s.init=s,z.init=z,K2D=K2D,J=J,K=K,X=X.new,n.marked=n.marked,
              xlim=xlim,ylim=ylim,locs=locs.use2,tel.inds=tel.inds,n.locs.ind=n.locs.ind,n.tel.inds=n.tel.inds,
              res=res,cellArea=cellArea,x.vals=x.vals,xlim=xlim,ylim=ylim,
              y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
              n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data))

}