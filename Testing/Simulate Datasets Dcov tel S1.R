setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/")
n.datasets <- 100
n.cores <- 18
seeds <- round(runif(n.datasets,1,10000))
source("sim.SMR.Dcov.R")

#simulate some data

####Simulate some data####
n.marked <- 20 #for "premarked", this is the number marked. for "natural", this is the number of captured individuals that are identifiable
lam0 <- 0.5
sigma <- 0.5
K <- 10 #number of occasions
buff <- 2 #state space buffer
X <- expand.grid(3:11,3:11) #make a trapping array
#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
obstype <- "poisson"
tlocs <- 0 #number of telemetry locs/marked individual. For "premarked"

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- range(X[,1]) + c(-buff,buff)
ylim <- range(X[,2]) + c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X[,1] <- X[,1]-x.shift
X[,2] <- X[,2]-y.shift

res <- 0.25 #habitat grid resolution, length of 1 cell side
cellArea <- res^2 #area of one cell
x.vals <- seq(xlim[1]+res/2,xlim[2]-res/2,res) #x cell centroids
y.vals <- seq(ylim[1]+res/2,ylim[2]-res/2,res) #y cell centroids
dSS <- as.matrix(cbind(expand.grid(x.vals,y.vals)))
cells <- matrix(1:nrow(dSS),nrow=length(x.vals),ncol=length(y.vals))
n.cells <- nrow(dSS)
n.cells.x <- length(x.vals)
n.cells.y <- length(y.vals)

#simulate a D.cov, higher cov.pars for large scale cov
#change seed to get new D.cov. trial and error to create one with good trapping array coverage
set.seed(152)
library(geoR)
D.cov <- grf(n.cells,grid=dSS,cov.pars=c(500,500),messages=FALSE)[[2]] #takes a while, run time depends on n.cells. 3600 cells pretty fast
D.cov <- as.numeric(scale(D.cov)) #scale
par(mfrow=c(1,1),ask=FALSE)
image(data$x.vals,data$y.vals,matrix(data$D.cov,data$n.cells.x,data$n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X,pch=4)

#Additionally, maybe we want to exclude "non-habitat"
#just removing the corners here for simplicity
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(1,length(D.cov))
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]<2&dSS.tmp[,2]>12] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]<2] <- 0
InSS[dSS.tmp[,1]>12&dSS.tmp[,2]>12] <- 0

image(x.vals,y.vals,matrix(InSS,n.cells.x,n.cells.y),main="Habitat")

#Density covariates
D.beta0 <- -0.5 #data simulator uses intercept for marked + unmarked
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X,pch=4,cex=1,lwd=2)


setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing")
dir.create("Dcov_S1_datasets")
setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Dcov_S1_datasets")
library(snow)
library(doSNOW)
library(foreach)
cl.tmp  <-  makeCluster(rep("localhost",n.cores), type="SOCK")
registerDoSNOW(cl.tmp)
out <- foreach(rep=1:n.datasets) %dopar% {
  data <- sim.SMR.Dcov(D.beta0=D.beta0,D.beta1=D.beta1,res=res,
                       D.cov=D.cov,InSS=InSS,n.marked=n.marked,
                       theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                       lam0=lam0,sigma=sigma,K=K,X=X,xlim=xlim,ylim=ylim,tlocs=tlocs,
                       obstype=obstype)
  save(data,file=paste0("S1_data_",rep,".RData"))
}
