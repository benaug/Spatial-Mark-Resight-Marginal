#This is an SMR data simulator and MCMC sampler that handles all sample types
#1) marked, known ID
#2) marked, unknown ID
#3) unmarked, unknown ID
#4) unknown marked status, unknown ID

#It also includes a marking process for "generalized SMR" or "spatial capture-mark-resight"
#with interspersed marking and sighting

#y[i,j,k] ~ Poisson(lam[i,j,k])
#y.event[i,j,k,1:3] ~ Multinomial(theta.marked[1:3],y[i,j,k]) for marked i
#y.event[i,j,k,1:3] ~ Multinomial(theta.unmarked[1:3],y[i,j,k]) for unmarked i

#event 1 is you know the ID (marked known ID samples)
#event 2 is you know the mark status, but not ID (marked, unknown ID or unmarked samples)
#event 3 is you don't know mark status or ID (unknown marked status samples)

library(nimble)
source("sim.SMR.Dcov.Generalized.Interspersed.R")
source("NimbleModel SMR Poisson Dcov DA2 Marginal Generalized Interspersed.R")
source("NimbleFunctions SMR Poisson Dcov DA2 Marginal Generalized Interspersed.R")
source("init.SMR.Dcov.Generalized.Interspersed.R")
source("sSampler Poisson Dcov Marginal Generalized Interspersed.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

####Simulate some data####
p0 <- 0.2 #marking process p0
lam0 <- 0.5 #sighting process lam0
sigma <- 0.5 #shared sigma
K.mark <- 5 #number of marking occasions
K.sight <- 10 #number of sighting occasions
buff <- 2 #state space buffer
X.mark <- expand.grid(5:9,5:9) #marking process traps
X.sight <- expand.grid(3:11,3:11) #sighting process traps
#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- c(0.75,0.15,0.1) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- 0.75 #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
obstype <- "poisson" #for sighting process. can also simulate "negbin", but cannot be fit with marginalized approach
tlocs <- 10 #number of telemetry locs/marked individual.

# a character vector of length K.mark+K.sight specifying the order of the marking and sighting occasions.
#Vector elements are either "M" or "S" arranged in the order the marking and sighting occasions occurred.
#there must be K.mark M's and K.sight S's. Data simulator will check these requirements.
K.order <- c("S","M","S","M","S","M","S","M","S","M","S","S","S","S","S")

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

### Habitat Covariate stuff###
#get x and y extent by buffering state space
X.both <- rbind(X.mark,X.sight)
xlim <- range(X.both[,1]) + c(-buff,buff)
ylim <- range(X.both[,2]) + c(-buff,buff)
#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
x.shift <- xlim[1]
y.shift <- ylim[1]
xlim <- xlim-x.shift
ylim <- ylim-y.shift
X.both[,1] <- X.both[,1]-x.shift
X.both[,2] <- X.both[,2]-y.shift
X.mark[,1] <- X.mark[,1]-x.shift
X.mark[,2] <- X.mark[,2]-y.shift
X.sight[,1] <- X.sight[,1]-x.shift
X.sight[,2] <- X.sight[,2]-y.shift


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
points(X.sight,pch=4,lwd=2)
points(X.mark,pch=4,col="darkred",lwd=2)

#Additionally, maybe we want to exclude "non-habitat" or limit the state space extent
#let's use a 3sigma buffer
dSS.tmp <- dSS - res/2 #convert back to grid locs
InSS <- rep(0,length(D.cov))
dists <- e2dist(X.both,dSS.tmp)
min.dists <- apply(dists,2,min)
InSS[min.dists<(3*sigma)] <- 1
image(x.vals,y.vals,matrix(D.cov*InSS,n.cells.x,n.cells.y),main="Habitat",col=cols1)
points(X.sight,pch=4,lwd=2)
points(X.mark,pch=4,col="darkred",lwd=2)


#Density covariates
D.beta0 <- -0.5 #data simulator uses intercept for marked + unmarked
D.beta1 <- 0.5
#what is implied expected N in state space?
lambda.cell <- InSS*exp(D.beta0 + D.beta1*D.cov)*cellArea
sum(lambda.cell) #expected N in state space

image(x.vals,y.vals,matrix(lambda.cell,n.cells.x,n.cells.y),main="Expected Density",col=cols1)
points(X.sight,pch=4,cex=1,lwd=2)
points(X.mark,pch=4,cex=1,lwd=2,col="darkred")

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(143532) #change seed for new data set
data <- sim.SMR.Dcov.Generalized.Interspersed(D.beta0=D.beta0,D.beta1=D.beta1,res=res,
                D.cov=D.cov,InSS=InSS,
                theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                lam0=lam0,p0=p0,sigma=sigma,
                K.mark=K.mark,K.sight=K.sight,K.order=K.order,
                X.mark=X.mark,X.sight=X.sight,xlim=xlim,ylim=ylim,tlocs=tlocs,
                obstype=obstype)
points(data$s,pch=16) #add activity centers
points(data$s.marked,col="goldenrod",pch=16) #marked individual activity centers
data$n.marked #number of individuals captured and marked

#What is the observed data?
str(data$y.mark) #marking process captures
str(data$y.mID) #marked with ID detections
str(data$y.mnoID) #marked with no ID samples
str(data$y.um) #unmarked samples
str(data$y.unk) #unknown marked status samples
str(data$marked.status) #history of which sighting occasions marked individuals are marked
str(data$locs) #possibly telemetry. n.marked x tlocs x 2 array (or ragged array if number of locs/ind differ). 
#Rows are 1:n.marked individuals, columns are max telemetry points for a single
#individual, fill in NAs for inds with no telemetry and/or inds without max number of telemetry points.
#in latter case, order telemetry points first, then NAs

#function to test for errors in mask set up. 
mask.check(dSS=data$dSS,cells=data$cells,n.cells=data$n.cells,n.cells.x=data$n.cells.x,
           n.cells.y=data$n.cells.y,res=data$res,xlim=data$xlim,ylim=data$ylim,
           x.vals=data$x.vals,y.vals=data$y.vals)

####Fit model in Nimble####
M <- 150
X.mark <- data$X.mark
J.mark <- nrow(X.mark)
K.mark <- data$K.mark
n.marked <- data$n.marked
X.sight <- data$X.sight
J.sight <- nrow(X.sight)
K.sight <- data$K.sight

#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#provide inits for theta.marked and theta.unmarked to initialize data. Assuming all sessions the same
#each must sum to 1
theta.marked.init <- c(0.5,0.25,0.25)
theta.unmarked.init <- c(0,0.75,0.25) #first cell for theta.unmarked.init should be 0
inits <- list(p0=0.5,lam0=1,sigma=1,theta.marked=theta.marked.init,theta.unmarked=theta.unmarked.init) #initializing with 1 parameter per session, just set all to same value


#This function structures the simulated data to fit the model in Nimble (some more restructing below)
#Also checks some inits
nimbuild <- init.SMR.Dcov.Generalized.Interspersed(data,inits,M=M)
image(data$x.vals,data$y.vals,matrix(data$D.cov*data$InSS,data$n.cells.x,data$n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
points(nimbuild$s,pch=16) #initialized activity centers
points(X.sight,pch=4,col="lightblue")
points(X.mark,pch=4)
y.mark2D <- apply(data$y.mark,c(1,2),sum)
y.mID2D <-  apply(data$y.mID,c(1,2),sum)
for(i in 1:n.marked){
  trapcaps1 <- which(y.mark2D[i,]>0)
  trapcaps2 <- which(data$y.mID2D[i,]>0)
  traps <-  rbind(X.mark[trapcaps1,],X.sight[trapcaps2,])
  s <- nimbuild$s[i,]
  points(s[1],s[2],col="goldenrod",pch=16)
  if(nrow(traps)>0){
    for(j in 1:nrow(traps)){
      lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
    }
  }
}

#inits for nimble
D0.init <- (sum(nimbuild$z))/(sum(data$InSS)*data$res^2)

#must initialize N.M and N.UM to be consistent with z. speeds converge to set consistent with lambda.N.M/UM
Niminits <- list(z=nimbuild$z,s=nimbuild$s,D0=D0.init,D.beta1=0,
                 N=sum(nimbuild$z),
                 theta.unmarked=c(0,0.5,0.5),
                 p0=inits$p0,lam0=inits$lam0,sigma=inits$sigma)

dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are

#Use this if you do not have telemetry. Make sure telemetry commented out in model file
# constants <- list(n.marked=n.marked,M=M,
#                   J.mark=J.mark,K1D.mark=data$K1D.mark,
#                   J.sight=J.sight,K2D.sight=data$K2D.sight,
#                   D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
#                   xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res)
# 
# Nimdata <- list(y.mark=nimbuild$y.mark,
#                 y.mID=nimbuild$y.mID, #marked with ID
#                 y.mnoID=nimbuild$y.mnoID, #marked without ID
#                 y.um=nimbuild$y.um, #unmarked
#                 y.unk=nimbuild$y.unk, #unk marked status
#                 marked.status=data$marked.status, #individual by occasion marked status for marked guys
#                 X.mark=as.matrix(X.mark),X.sight=as.matrix(X.sight),
#                 dummy.data=dummy.data,cells=data$cells,InSS=data$InSS)

#If you have telemetry use these instead. Make sure to uncomment telemetry BUGS code.
constants <- list(n.marked=n.marked,M=M,
                  J.mark=J.mark,K1D.mark=data$K1D.mark,
                  J.sight=J.sight,K.sight=data$K.sight,K2D.sight=data$K2D.sight,
                  D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res,
                  tel.inds=nimbuild$tel.inds,
                  n.tel.inds=length(nimbuild$tel.inds),n.locs.ind=nimbuild$n.locs.ind)
Nimdata <- list(y.mark=nimbuild$y.mark, #marking process
                y.mID=nimbuild$y.mID, #marked with ID
                y.mnoID=nimbuild$y.mnoID, #marked without ID
                y.um=nimbuild$y.um, #unmarked
                y.unk=nimbuild$y.unk, #unk marked status
                marked.status=data$marked.status, #individual by occasion marked status for marked guys
                dummy.data=dummy.data,cells=data$cells,InSS=data$InSS,
                X.mark=as.matrix(data$X.mark),X.sight=as.matrix(data$X.sight),locs=data$locs)

# set parameters to monitor
parameters <- c('D0','lambda.N','p0','lam0','sigma','theta.marked','theta.unmarked',
                'N','D.beta1')
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#use block sampler below for 'D0.M','D0.UM','D.beta1'
config.nodes <- c('p0','lam0','sigma','theta.marked','theta.unmarked[2:3]')
# config.nodes <- c()
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                      monitors2=parameters2,thin2=10,nodes=config.nodes)

# how many z proposals per iteration per session for marked (if updated), unmarked?
z.ups <- round(M*0.25) #doing 25% of M
#nodes used for update, calcNodes + z nodes
pd.nodes <- Rmodel$expandNodeNames("pd")
lam.nodes <- Rmodel$expandNodeNames("lam")
y.mark.nodes <- Rmodel$expandNodeNames("y.mark")
y.um.nodes <- Rmodel$expandNodeNames("y.um")
y.unk.nodes <- Rmodel$expandNodeNames("y.unk")
bigLam.unmarked.nodes <- Rmodel$expandNodeNames("bigLam.unmarked") #only need this in calcNodes
lam.um.nodes <- Rmodel$expandNodeNames("lam.um")
lam.unk.nodes <- Rmodel$expandNodeNames("lam.unk")
N.node <- Rmodel$expandNodeNames("N")
z.nodes <- Rmodel$expandNodeNames("z")
calcNodes <- c(N.node,pd.nodes,lam.nodes,bigLam.unmarked.nodes,
               lam.um.nodes,lam.unk.nodes,y.mark.nodes,y.um.nodes,y.unk.nodes)
conf$addSampler(target = paste("N"),
                type = 'zSampler',control = list(z.ups=z.ups,J.sight=J.sight,K.sight=K.sight,
                                                 n.marked=n.marked,M=M,N.node=N.node,
                                                 pd.nodes=pd.nodes,lam.nodes=lam.nodes,
                                                 lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                 y.mark.nodes=y.mark.nodes,y.um.nodes=y.um.nodes,
                                                 y.unk.nodes=y.unk.nodes,
                                                 z.nodes=z.nodes,calcNodes=calcNodes),
                silent = TRUE)

#add sSampler
# if no telemetry,
# for(i in 1:M){
#   conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
#                   type = 'sSampler',control=list(i=i,J.mark=J.mark,J.sight=J.sight,K.sight=K.sight,res=data$res,
#                                                  n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
#                                                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,
#                                                  n.marked=n.marked,n.locs.ind=0,
#                                                  scale=0.25),silent = TRUE)
#   #scale parameter here is just the starting scale. It will be tuned.
# }
#with telemetry (make sure you turn it on in model code),
for(i in 1:M){
  if(i %in% nimbuild$tel.inds){#inds with telemetry
    conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,J.mark=J.mark,J.sight=J.sight,K.sight=K.sight,res=data$res,
                                                   n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                   xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                                                   n.locs.ind=nimbuild$n.locs.ind[i],
                                                   n.marked=n.marked,scale=0.25),silent = TRUE)
  }else{ #inds with no telemetry
    conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,J.mark=J.mark,J.sight=J.sight,K.sight=K.sight,res=data$res,
                                                   n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                   xlim=nimbuild$xlim,ylim=nimbuild$ylim,n.locs.ind=0,
                                                   n.marked=n.marked,scale=0.25),silent = TRUE)
  }
}

#can add block sampler if lam0, sigma, and/or lambda.N posteriors highly correlated
#probably not needed with telemetry
conf$addSampler(target = c("lam0","sigma"),type = 'RW_block',
                control = list(adaptive=TRUE),silent = TRUE)
#AF_slice pretty fast here
conf$addSampler(target = c("D0","D.beta1"),
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
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples <- as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[250:nrow(mvSamples),]))

exp(D.beta0)
data$N

tmp <- cor(mvSamples[250:nrow(mvSamples),])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)

#Important! If N hits M during sampling, raise M. 

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
