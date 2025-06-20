setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/")
n.datasets <- 100
n.cores <- 18
seeds <- round(runif(n.datasets,1,10000))
source("sim.SMR.multisession.Dcov.R")
source("sim.SMR.Dcov.R") #used by multisession simulator

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

####Simulate some data####
#Here, I'll simulate 3 populations with different n.marked, K, X, and state space areas
#sharing lam0, sigma so they can be shared during estimation
N.session <- 3
n.marked <- c(20,21,22)
lam0 <- rep(0.5,N.session)
sigma <- rep(0.5,N.session)
K <- c(10,11,12) #number of occasions
buff <- rep(2,N.session) #state space buffer

#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- matrix(rep(c(0.75,0.15,0.1),N.session),nrow=N.session,byrow=TRUE) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- rep(0.75,N.session) #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
obstype <- "poisson"
tlocs <- c(25,25,25) #number of telemetry locs/marked individual in each session. For "premarked"

#make an SCR trapping array. using same trapping array and Dcov here for each session
X <- vector("list",N.session)
X[[1]] <- expand.grid(3:11,3:11)
X[[2]] <- expand.grid(3:11,3:11)
X[[3]] <- expand.grid(3:11,3:11)

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- ylim <- matrix(NA,N.session,2)
for(g in 1:N.session){
  xlim[g,] <- range(X[[g]][,1]) + c(-buff[g],buff[g])
  ylim[g,] <- range(X[[g]][,2]) + c(-buff[g],buff[g])
}

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
for(g in 1:N.session){
  x.shift <- xlim[g,1]
  y.shift <- ylim[g,1]
  xlim[g,] <- xlim[g,] - x.shift
  ylim[g,] <- ylim[g,] - y.shift
  X[[g]][,1] <- X[[g]][,1]- x.shift
  X[[g]][,2] <- X[[g]][,2]- y.shift
}

res <- rep(0.25,N.session) #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- y.vals <- dSS <- cells <- vector("list",N.session)
n.cells <- n.cells.x <- n.cells.y <- rep(NA,N.session)
for(g in 1:N.session){
  x.vals[[g]] <- seq(xlim[g,1]+res[g]/2,xlim[g,2]-res[g]/2,res[g]) #x cell centroids
  y.vals[[g]] <- seq(ylim[g,1]+res[g]/2,ylim[g,2]-res[g]/2,res[g]) #y cell centroids
  dSS[[g]] <- as.matrix(cbind(expand.grid(x.vals[[g]],y.vals[[g]])))
  cells[[g]] <- matrix(1:nrow(dSS[[g]]),nrow=length(x.vals[[g]]),ncol=length(y.vals[[g]]))
  n.cells[g] <- nrow(dSS[[g]])
  n.cells.x[g] <- length(x.vals[[g]])
  n.cells.y[g] <- length(y.vals[[g]])
}

#create a density covariate - one for each session
library(geoR)
D.cov <- vector("list",N.session)
#need a simulated landscape with individuals living around traps to be captured
#these are pretty good
# D.seeds <- c(13223,13216,13252)
D.seeds <- c(152,152,152)
for(g in 1:N.session){
  set.seed(D.seeds[g])
  D.cov.tmp <- grf(n.cells[g],grid=dSS[[g]],cov.pars=c(500,500),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
  D.cov.tmp <- as.numeric(scale(D.cov.tmp)) #scale
  par(mfrow=c(1,1),ask=FALSE)
  D.cov[[g]] <- D.cov.tmp
  image(x.vals[[g]],y.vals[[g]],matrix(D.cov[[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g," D.cov"),xlab="X",ylab="Y",col=cols1)
}

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
InSS <- vector("list",N.session)
for(g in 1:N.session){
  dSS.tmp <- dSS[[g]] - res[g]/2 #convert back to grid locs
  InSS[[g]] <- rep(1,length(D.cov[[g]]))
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]<2&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]<2] <- 0
  InSS[[g]][dSS.tmp[,1]>(xlim[g,2]-2)&dSS.tmp[,2]>(ylim[g,2]-2)] <- 0
  image(x.vals[[g]],y.vals[[g]],matrix(InSS[[g]],n.cells.x[g],n.cells.y[g]),main=paste("Session",g," Habitat"))
}

#Density covariates
D.beta0 <- rep(-0.5,N.session)
D.beta1 <- rep(0.5,N.session)
#what is implied expected N in state space?
for(g in 1:N.session){
  lambda.cell <- InSS[[g]]*exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  print(sum(lambda.cell)) #expected N in state space
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell,n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g," Expected Density"),col=cols1)
  points(X[[g]],pch=4,cex=0.75)
}

setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing")
dir.create("Multisession_tel_S1_datasets")
setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Multisession_tel_S1_datasets")
library(snow)
library(doSNOW)
library(foreach)
cl.tmp  <-  makeCluster(rep("localhost",n.cores), type="SOCK")
registerDoSNOW(cl.tmp)
out <- foreach(rep=1:n.datasets) %dopar% {
  data <- sim.SMR.multisession.Dcov(N.session=N.session,n.marked=n.marked,
                                    theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                                    lam0=lam0,sigma=sigma,K=K,X=X,tlocs=tlocs,obstype=obstype,
                                    D.beta0=D.beta0,D.beta1=D.beta1,D.cov=D.cov,InSS=InSS,res=res,xlim=xlim,ylim=ylim)
  save(data,file=paste0("S1_data_",rep,".RData"))
}

