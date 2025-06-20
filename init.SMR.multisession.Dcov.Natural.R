e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

init.SMR.multisession.Dcov.Natural <- function(data,inits=NA,M1=NA,M2=NA){
  N.session <- length(data)
  n.marked <- sapply(data,function(x){x$n.marked})
  if(length(M1)!=N.session)stop("Must supply an M1 for each session.")
  if(length(M2)!=N.session)stop("Must supply an M2 for each session.")
  init.session <- vector("list",N.session)
  for(g in 1:N.session){
    inits.use <- inits #lam0 and sigma inits vary by session
    inits.use$lam0 <- inits.use$lam0[g]
    inits.use$sigma <- inits.use$sigma[g]
    init.session[[g]] <- init.SMR.Dcov.Natural(data[[g]],inits.use,M1=M1[g],M2=M2[g])
  }
  
  J <- unlist(lapply(data,function(x){nrow(x$X)}))
  M.both <- M1 + M2
  maxM <- max(M.both)
  s <- array(NA,dim=c(N.session,maxM,2))
  z <- matrix(NA,N.session,maxM)
  K1D <- matrix(NA,N.session,max(J))
  y.mID <- array(NA,dim=c(N.session,max(M1),max(J)))
  y.mnoID <- matrix(NA,N.session,max(J))
  y.um <- matrix(NA,N.session,max(J))
  y.unk <- matrix(NA,N.session,max(J))
  
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
    s[g,1:M.both[g],] <- init.session[[g]]$s
    z[g,1:M.both[g]] <- init.session[[g]]$z
    K1D[g,1:J[g]] <- init.session[[g]]$K1D
    y.mID[g,1:M1[g],1:J[g]] <- init.session[[g]]$y.mID
    y.mnoID[g,1:J[g]] <- init.session[[g]]$y.mnoID
    y.um[g,1:J[g]] <- init.session[[g]]$y.um
    y.unk[g,1:J[g]] <- init.session[[g]]$y.unk
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
  X.new <- array(NA,dim=c(N.session,max(J),2))
  for(g in 1:N.session){
    X.new[g,1:J[g],] <- data[[g]]$X
  }
  dummy.data <- matrix(0,N.session,maxM) #dummy data not used, doesn't really matter what the values are
  
  
  return(list(y.mID=y.mID,y.mnoID=y.mnoID,y.um=y.um,y.unk=y.unk,
              s.init=s,z.init=z,K1D=K1D,J=J,X=X.new,n.marked=n.marked,
              xlim=xlim,ylim=ylim,
              res=res,cellArea=cellArea,x.vals=x.vals,xlim=xlim,ylim=ylim,
              y.vals=y.vals,dSS=dSS,InSS=InSS,cells=cells,n.cells=n.cells,n.cells.x=n.cells.x,
              n.cells.y=n.cells.y,D.cov=D.cov,dummy.data=dummy.data))
  
}