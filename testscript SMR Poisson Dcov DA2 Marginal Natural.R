#This is an SMR data simulator and MCMC sampler that handles all sample types
#1) marked, known ID
#2) marked, unknown ID
#3) unmarked, unknown ID
#4) unknown marked status, unknown ID

#It handles natural marks scenarios. marked, unknown ID detections allowed here, but
#these may or may not occur in practice. e.g., if you define "marked" as "has scar"
#and all scars can provide individual IDs, you can see a scar but not be able to make
#the ID. This makes sense. But if "marked" means it is an individual you have identified,
#it makes less sense that you can tell an identifiable individual is in a photo, particularly
#if you never obtained a marked with ID sample from this individual. Think about whether these should be
#unknown marked status samples.

#y[i,j,k] ~ Poisson(lam[i,j,k])
#y.event[i,j,k,1:3] ~ Multinomial(theta.marked[1:3],y[i,j,k]) for marked i
#y.event[i,j,k,1:3] ~ Multinomial(theta.unmarked[1:3],y[i,j,k]) for unmarked i

#event 1 is you know the ID (marked known ID samples)
#event 2 is you know the mark status, but not ID (marked, unknown ID or unmarked samples)
#event 3 is you don't know mark status or ID (unknown marked status samples)

library(nimble)
source("sim.SMR.Dcov.Natural.R")
source("NimbleModel SMR Poisson Dcov DA2 Marginal Natural.R")
source("NimbleFunctions SMR Poisson Dcov DA2 Marginal Natural.R")
source("init.SMR.Dcov.Natural.R")
source("sSampler Poisson Dcov Marginal Natural.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
lam0 <- 0.5
sigma <- 0.5
K <- 10 #number of occasions
buff <- 2 #state space buffer
X <- expand.grid(3:11,3:11) #make a trapping array
#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
obstype <- "poisson" #can also simulate "negbin", but cannot fit with marginalized approach
p.marked <- 0.25 #probability a detected individual is "marked"

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
image(x.vals,y.vals,matrix(D.cov,n.cells.x,n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X,pch=4)

#Additionally, maybe we want to exclude "non-habitat" or limit the state space extent
#let's use a 3sigma buffer
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(0,length(D.cov))
dists <- e2dist(X,dSS.tmp)
min.dists <- apply(dists,2,min)
InSS[min.dists<(3*sigma)] <- 1
image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="Habitat",col=cols1)
points(X,pch=4,col="darkred",lwd=2)

#Density covariates
D.beta0 <- -0.75 #data simulator uses intercept for marked + unmarked
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X,pch=4,cex=1,lwd=2)

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(143532) #change seed for new data set
data <- sim.SMR.Dcov.Natural(D.beta0=D.beta0,D.beta1=D.beta1,res=res,
                D.cov=D.cov,InSS=InSS,p.marked=p.marked,
                theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                lam0=lam0,sigma=sigma,K=K,X=X,xlim=xlim,ylim=ylim,
                obstype=obstype)
points(data$s,pch=16) #add activity centers


#What is the observed data?
str(data$y.mID) #marked with ID detections
str(data$y.mnoID) #marked with no ID samples
str(data$y.um) #unmarked samples
str(data$y.unk) #unknown marked status samples
str(data$locs) #possibly telemetry. n.marked x tlocs x 2 array (or ragged array if number of locs/ind differ). 
#Rows are 1:n.marked individuals, columns are max telemetry points for a single
#individual, fill in NAs for inds with no telemetry and/or inds without max number of telemetry points.
#in latter case, order telemetry points first, then NAs

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

####Fit model in Nimble####
n.marked <- data$n.marked
n.marked #M1 needs to be larger than number of marked inds detected
M1 <- 75
M2 <- 100
M.both <- M1 + M2
X <- data$X
J <- nrow(X)
K <- data$K

#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#also use this function checks to make sure theta.marked and theta.unmarked inits are in
#the correct structure. 
inits <- list(lam0=1,sigma=1)

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.SMR.Dcov.Natural(data,inits,M1=M1,M2=M2)
image(data$x.vals,data$y.vals,matrix(data$D.cov,data$n.cells.x,data$n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(X,pch=4)
points(nimbuild$s,pch=16) #initialized activity centers

#inits for nimble
D0.M.init <- (sum(nimbuild$z[1:M1]))/(sum(data$InSS)*data$res^2)
D0.UM.init <- (sum(nimbuild$z[(M1+1):M.both]))/(sum(data$InSS)*data$res^2)

#must initialize N.M and N.UM to be consistent with z. speeds converge to set consistent with lambda.N.M/UM
Niminits <- list(z=nimbuild$z,s=nimbuild$s,D0.M=D0.M.init,D0.UM=D0.UM.init,D.beta1=0,
                 N.M=sum(nimbuild$z[1:M1]),N.UM=sum(nimbuild$z[(M1+1):M.both]),
                 theta.unmarked=c(0,0.5,0.5),lam0=inits$lam0,sigma=inits$sigma)

dummy.data <- rep(0,M.both) #dummy data not used, doesn't really matter what the values are

#Use this if you do not have telemetry. Make sure telemetry commented out in model file
constants <- list(M1=M1,M.both=M.both,J=J,K1D=data$K1D,
                  D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res)

Nimdata <- list(y.mID=nimbuild$y.mID, #marked with ID
                y.mnoID=nimbuild$y.mnoID, #marked without ID
                y.um=nimbuild$y.um, #unmarked
                y.unk=nimbuild$y.unk, #unk marked status
                X=as.matrix(X),
                dummy.data=dummy.data,cells=data$cells,InSS=data$InSS)

# set parameters to monitor
parameters <- c('D0','D.beta1','lambda.N.M','lambda.N.UM','N','N.UM',"N.M",
                'lam0','sigma','theta.marked','theta.unmarked')
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#use block sampler below for 'D0.M','D0.UM','D.beta1'
config.nodes <- c('lam0','sigma','theta.marked','theta.unmarked[2:3]')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                      monitors2=parameters2,thin2=10,nodes=config.nodes)

# how many z proposals per iteration per session for marked (if updated), unmarked?
z.ups <- round(c(M1*0.25,M2*0.25)) #doing 25% of M1 and M2 here
z.ups #number of N/z ups for marked and unmarked inds per iteration
#nodes used for update, calcNodes + z nodes
lam.nodes <- Rmodel$expandNodeNames(paste("lam[1:",M.both,",1:",J,"]"))
y.mID.nodes <- Rmodel$expandNodeNames(paste("y.mID[1:",M1,",1:",J,"]"))
y.mnoID.nodes <- Rmodel$expandNodeNames(paste("y.mnoID[1:",J,"]"))
y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[1:",J,"]"))
y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[1:",J,"]"))
bigLam.marked.nodes <- Rmodel$expandNodeNames("bigLam.marked") #only need this in calcNodes
bigLam.unmarked.nodes <- Rmodel$expandNodeNames("bigLam.unmarked") #only need this in calcNodes
lam.mnoID.nodes <- Rmodel$expandNodeNames("lam.mnoID")
lam.um.nodes <- Rmodel$expandNodeNames("lam.um")
lam.unk.nodes <- Rmodel$expandNodeNames("lam.unk")
N.M.node <- Rmodel$expandNodeNames(paste("N.M"))
N.UM.node <- Rmodel$expandNodeNames(paste("N.UM"))
N.node <- Rmodel$expandNodeNames(paste("N"))
z.nodes <- Rmodel$expandNodeNames(paste("z[1:",M.both,"]"))
calcNodes <- c(N.M.node,N.UM.node,N.node=N.node,lam.nodes,bigLam.marked.nodes,bigLam.unmarked.nodes,
               lam.mnoID.nodes,lam.um.nodes,lam.unk.nodes,z.nodes,
               y.mID.nodes,y.mnoID.nodes,y.um.nodes,y.unk.nodes)
conf$addSampler(target = paste("N.UM"),
                type = 'zSampler',control = list(z.ups=z.ups,J=J,n.marked=n.marked,M1=M1,M.both=M.both,
                                                 N.M.node=N.M.node,N.UM.node=N.UM.node,
                                                 lam.nodes=lam.nodes,lam.mnoID.nodes=lam.mnoID.nodes,
                                                 lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                 y.mID.nodes=y.mID.nodes,y.mnoID.nodes=y.mnoID.nodes,
                                                 y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                                 z.nodes=z.nodes,calcNodes=calcNodes),
                silent = TRUE)

#add sSampler
for(i in 1:M.both){
  conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                  type = 'sSampler',control=list(i=i,J=J,M1=M1,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                 xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                                                 scale=0.25),silent = TRUE)
  #scale parameter here is just the starting scale. It will be tuned.
}

#can add block sampler if lam0, sigma, and/or lambda.N posteriors highly correlated
conf$addSampler(target = c("lam0","sigma"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
# AF_slice pretty fast here
conf$addSampler(target = c("D0.M","D0.UM","D.beta1"),
                type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=1) #this will run in R, used for better debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# Run the model.
start.time2 <- Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time <- Sys.time()
end.time - start.time  # total time for compilation, replacing samplers, and fitting
end.time - start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[200:nrow(mvSamples),]))

exp(D.beta0)
data$N

tmp <- cor(mvSamples[250:nrow(mvSamples),])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)

#Important! If N.M hits M1 during sampling, raise M1. If N.UM hits M2 during sampling, raise M2.

#plot density surface, etc.
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10

#compare expected D plot to truth (for simulated data sets)
n.cells <- data$n.cells
lambda.cell <- exp(D.beta0 + D.beta1*D.cov)*cellArea
n.iter.use <- burnin2:nrow(mvSamples2)
lambda.cell.post <- t(cellArea*mvSamples2[n.iter.use,D0.idx]*mvSamples2[n.iter.use,lambda.cell.idx[1:n.cells]])
lambda.cell.ests <- rowMeans(lambda.cell.post[1:n.cells,])
lambda.cell.HPDs <- HPDinterval(mcmc(t(lambda.cell.post[1:n.cells,])))
#remove nonhabitat (or not, comment out)
lambda.cell[data$InSS==0] <- NA
lambda.cell.ests[data$InSS==0] <- NA

par(mfrow=c(1,1),ask=FALSE)
zlim <- range(c(lambda.cell,lambda.cell.ests),na.rm=TRUE) #use same zlim for plots below
#truth
image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)
#estimate, posterior means
image(x.vals,y.vals,matrix(lambda.cell.ests,n.cells.x,n.cells.y),main="Expected Density",zlim=zlim)

#cell ests and 95% HPDs vs. truth. 
#Need a lot of posterior samples for accurate 95% HPDs, if not, will look "jagged"
idx <- order(lambda.cell)
plot(lambda.cell.ests[1:n.cells][idx]~lambda.cell[1:n.cells][idx],type="l",lwd=2,
     main="True vs. Estimated Density",ylim=range(lambda.cell.HPDs[1:n.cells,]))
lines(lambda.cell.HPDs[1:n.cells,1][idx]~lambda.cell[1:n.cells][idx],lty=2)
lines(lambda.cell.HPDs[1:n.cells,2][idx]~lambda.cell[1:n.cells][idx],lty=2)
abline(0,1,col="darkred",lwd=2) #1:1 expectation

