e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.multisession.Dcov.Generalized <- function(data,inits=NA,M=NA){
  N.session <- length(data)
  n.marked <- sapply(data,function(x){x$n.marked})
  if(length(M)!=N.session)stop("Must supply an M for each session.")
  init.session <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #p0, lam0, and sigma inits vary by session
    inits.use$p0 <- inits.use$p0[g]
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.session[[g]] <- init.SMR.Dcov.Generalized(data[[g]],inits.use,M=M[g])
  }
  anyTelemetry <- FALSE
  n.tel.inds <- 0
  for(g in 1:N.session){
    if(!all(is.na(init.session[[g]]$tel.ID))){
      n.tel.inds <- n.tel.inds + length(init.session[[g]]$tel.ID)
      anyTelemetry <- TRUE
    }
  }
  if(anyTelemetry){
    tel.ID <- rep(NA,n.tel.inds) #population ID within session for each telemetry row
    tel.session <- rep(NA,n.tel.inds) #session for each telemetry row
    n.locs.ind <- rep(NA,n.tel.inds) #number of locations for each telemetry row
    max.tel.locs <- 0
    for(g in 1:N.session){
      if(!all(is.na(init.session[[g]]$tel.ID))){
        this.max <- max(init.session[[g]]$n.locs.ind)
        if(this.max>max.tel.locs){
          max.tel.locs <- this.max
        }
      }
    }
    locs.use2 <- array(NA,dim=c(n.tel.inds,max.tel.locs,2))
    ii <- 0
    for(g in 1:N.session){
      if(!all(is.na(init.session[[g]]$tel.ID))){
        for(t in 1:length(init.session[[g]]$tel.ID)){
          ii <- ii + 1
          tel.ID[ii] <- init.session[[g]]$tel.ID[t] 
          tel.session[ii] <- g 
          n.locs.ind[ii] <- init.session[[g]]$n.locs.ind[t] 
          locs.use2[ii,1:n.locs.ind[ii],] <- init.session[[g]]$locs[t,1:n.locs.ind[ii],]
        }
      }
    }
  }else{
    tel.ID <- tel.session <- n.locs.ind <- n.tel.inds <- NA
    locs.use2 <- NA
  }
  
  J.mark <- unlist(lapply(data,function(x){nrow(x$X.mark)}))
  J.sight <- unlist(lapply(data,function(x){nrow(x$X.sight)}))
  n.marked <- unlist(lapply(data,function(x){x$n.marked}))
  maxM <- max(M)
  s <- array(NA,dim=c(N.session,maxM,2))
  z <- matrix(NA,N.session,maxM)
  K1D.mark <- matrix(NA,N.session,max(J.mark))
  K1D.sight <- matrix(NA,N.session,max(J.sight))
  y.mark <- array(NA,dim=c(N.session,maxM,max(J.mark)))
  y.mID <- array(NA,dim=c(N.session,maxM,max(J.sight)))
  y.mnoID <- matrix(NA,N.session,max(J.sight))
  y.um <- matrix(NA,N.session,max(J.sight))
  y.unk <- matrix(NA,N.session,max(J.sight))
  
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
    K1D.mark[g,1:J.mark[g]] <- init.session[[g]]$K1D.mark
    K1D.sight[g,1:J.sight[g]] <- init.session[[g]]$K1D.sight
    y.mark[g,1:M[g],1:J.mark[g]] <- init.session[[g]]$y.mark
    y.mID[g,1:n.marked[g],1:J.sight[g]] <- init.session[[g]]$y.mID
    y.mnoID[g,1:J.sight[g]] <- init.session[[g]]$y.mnoID
    y.um[g,1:J.sight[g]] <- init.session[[g]]$y.um
    y.unk[g,1:J.sight[g]] <- init.session[[g]]$y.unk
    xlim[g,] <- data[[g]]$xlim
    ylim[g,] <- data[[g]]$ylim
    x.vals[g,1:n.cells.x[g]] <- data[[g]]$x.vals
    y.vals[g,1:n.cells.y[g]] <- data[[g]]$y.vals
    dSS[g,1:n.cells[g],] <- data[[g]]$dSS
    InSS[g,1:n.cells[g]] <- data[[g]]$InSS
    D.cov[g,1:n.cells[g]] <- data[[g]]$D.cov
    cells[g,1:n.cells.x[g],1:n.cells.y[g]] <- data[[g]]$cells
  }
  
  #put X in ragged array
  X.mark.new <- array(NA,dim=c(N.session,max(J.mark),2))
  X.sight.new <- array(NA,dim=c(N.session,max(J.sight),2))
  for(g in 1:N.session){
    X.mark.new[g,1:J.mark[g],] <- data[[g]]$X.mark
    X.sight.new[g,1:J.sight[g],] <- data[[g]]$X.sight
  }
  dummy.data <- matrix(0,N.session,maxM) #dummy data not used, doesn't really matter what the values are
  
  
  return(list(y.mark=y.mark,y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              s.init=s,z.init=z,K1D.mark=K1D.mark,K1D.sight=K1D.sight,J.mark=J.mark,J.sight=J.sight,
              X.mark=X.mark.new,X.sight=X.sight.new,n.marked=n.marked,
              xlim=xlim,ylim=ylim,locs=locs.use2,tel.ID=tel.ID,tel.session=tel.session,
              n.locs.ind=n.locs.ind,n.tel.inds=n.tel.inds, 
              res=res,cellArea=cellArea,x.vals=x.vals,xlim=xlim,ylim=ylim,
              y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
              n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data))
  
}