#multisession version of Poisson Dcov Marginal Generalized

library(nimble)
source("sim.SMR.multisession.Dcov.Generalized.R")
source("sim.SMR.Dcov.Generalized.R") #used by multisession simulator
source("NimbleModel SMR Multisession Poisson Dcov Marginal Generalized.R")
source("NimbleFunctions SMR Multisession Poisson Dcov Marginal Generalized.R")
source("init.SMR.multisession.Dcov.Generalized.R")
source("init.SMR.Dcov.Generalized.R") #used by multisession initializer
source("sSampler Multisession Poisson Dcov Marginal Generalized.R")
source("mask.check.R")

#If using Nimble version 0.13.1 and you must run this line 
nimbleOptions(determinePredictiveNodesInModel = FALSE)
# #If using Nimble before version 0.13.1, run this line instead
# nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

####Simulate some data####
N.session <- 3 #number of sessions
p0 <- rep(0.05,N.session) #marking process p0
lam0 <- rep(0.5,N.session) #sighting process lam0
sigma <- rep(0.5,N.session) #shared sigma
K.mark <- c(4,5,6) #number of marking occasions
K.sight <- c(10,9,8) #number of sighting occasions
buff <- rep(2,N.session) #state space buffer

#theta is probability of observing each sample type for marked and unmarked individuals
theta.marked <- matrix(rep(c(0.75,0.15,0.1),N.session),nrow=N.session,byrow=TRUE) #P(ID, Marked no ID, unk status). must sum to 1
theta.unmarked <- rep(0.75,N.session) #prob known marked status. #P(ID, Marked no ID, unk status)=(0,theta.unmarked,1-theta.unmarked)
obstype <- "poisson" #can also simulate "negbin", but cannot fit with marginalized approach
tlocs <- c(10,10,10) #number of telemetry locs/marked individual in each session.

X.mark <- vector("list",N.session)
X.mark[[1]] <- as.matrix(expand.grid(5:9,5:9))
X.mark[[2]] <- as.matrix(expand.grid(5:9,5:9))
X.mark[[3]] <- as.matrix(expand.grid(5:9,5:9))
X.sight <- vector("list",N.session)
X.sight[[1]] <- as.matrix(expand.grid(3:11,3:11))
X.sight[[2]] <- as.matrix(expand.grid(3:11,3:11))
X.sight[[3]] <- as.matrix(expand.grid(3:11,3:11))

#get some colors
library(RColorBrewer)
cols1 <- brewer.pal(9,"Greens")
cols2 <- brewer.pal(9,"YlOrBr")

### Habitat Covariate stuff###
#get x and y extent by buffering state space
xlim <- ylim <- matrix(NA,N.session,2)
X.both <- vector("list",N.session)
for(g in 1:N.session){
  X.both[[g]] <- rbind(X.mark[[g]],X.sight[[g]])
  xlim[g,] <- range(X.both[[g]][,1]) + c(-buff[g],buff[g])
  ylim[g,] <- range(X.both[[g]][,2]) + c(-buff[g],buff[g])
}

#shift X, xlim, ylim, so lower left side of state space is (0,0)
#this is required to use efficient look-up table to find the cell number
#of a continuous location
for(g in 1:N.session){
  x.shift <- xlim[g,1]
  y.shift <- ylim[g,1]
  xlim[g,] <- xlim[g,] - x.shift
  ylim[g,] <- ylim[g,] - y.shift
  X.both[[g]][,1] <- X.both[[g]][,1]- x.shift
  X.both[[g]][,2] <- X.both[[g]][,2]- y.shift
  X.mark[[g]][,1] <- X.mark[[g]][,1]- x.shift
  X.mark[[g]][,2] <- X.mark[[g]][,2]- y.shift
  X.sight[[g]][,1] <- X.sight[[g]][,1]- x.shift
  X.sight[[g]][,2] <- X.sight[[g]][,2]- y.shift
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
D.seeds <- c(13223,13216,13252)
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
  image(x.vals[[g]],y.vals[[g]],matrix(InSS[[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g," Habitat"),col=cols1)
  points(X.sight[[g]],pch=4,lwd=2)
  points(X.mark[[g]],pch=4,col="darkred",lwd=2)
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
  points(X.sight[[g]],pch=4,cex=1,lwd=2)
  points(X.mark[[g]],pch=4,cex=1,lwd=2,col="darkred")
}

#Simulate some data
#setting seed here because I am setting a seed to produce the D.cov and you will simulate the same
#data set over and over if you don't use different seeds here for each data set you simulate
set.seed(143532) #change seed for new data set
data <- sim.SMR.multisession.Dcov.Generalized(N.session=N.session,D.beta0=D.beta0,
                D.beta1=D.beta1,res=res,D.cov=D.cov,InSS=InSS,
                theta.marked=theta.marked,theta.unmarked=theta.unmarked,
                lam0=lam0,p0=p0,sigma=sigma,K.mark=K.mark,K.sight=K.sight,
                X.mark=X.mark,X.sight=X.sight,xlim=xlim,ylim=ylim,tlocs=tlocs,
                obstype=obstype)

#What is the observed data?
g <- 1 #session to look at
#What is the observed data?
str(data[[g]]$y.mark) #marking process captures
str(data[[g]]$y.mID) #marked with ID detections
str(data[[g]]$y.mnoID) #marked with no ID samples
str(data[[g]]$y.um) #unmarked samples
str(data[[g]]$y.unk) #unknown marked status samples
str(data[[g]]$locs) #possibly telemetry. n.marked x tlocs x 2 array (or ragged array if number of locs/ind differ). 

#Visualize activity centers
for(g in 1:N.session){
  lambda.cell <- InSS[[g]]*exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  image(data[[g]]$x.vals,data[[g]]$y.vals,matrix(lambda.cell,data[[g]]$n.cells.x,data[[g]]$n.cells.y),
        main=paste("Session",g," Expected Density"),col=cols1)
  points(X.sight[[g]],pch=4,cex=1,lwd=2)
  points(X.mark[[g]],pch=4,cex=1,lwd=2,col="darkred")
  points(data[[g]]$s,pch=16)
  points(data[[g]]$s.marked,col="goldenrod",pch=16) #marked individual activity centers
}

for(g in 1:N.session){
  #function to test for errors in mask set up. 
  mask.check(dSS=data[[g]]$dSS,cells=data[[g]]$cells,n.cells=data[[g]]$n.cells,n.cells.x=data[[g]]$n.cells.x,
             n.cells.y=data[[g]]$n.cells.y,res=data[[g]]$res,xlim=data[[g]]$xlim,ylim=data[[g]]$ylim,
             x.vals=data[[g]]$x.vals,y.vals=data[[g]]$y.vals)
}


####Fit model in Nimble####
M <- c(150,150,150)
X.mark <- lapply(data,function(x){x$X.mark})
X.sight <- lapply(data,function(x){x$X.sight})
J.mark <- sapply(X.mark,nrow)
J.sight <- sapply(X.sight,nrow)
K.mark <- sapply(data,function(x){x$K.mark})
K.sight <- sapply(data,function(x){x$K.sight})
n.marked <- sapply(data,function(x){x$n.marked})

#Need some inits to initialize data
#Use reasonable inits for lam0 and sigma since we check to make sure initial observation
#model likelihood is finite
#provide inits for theta.marked and theta.unmarked to initialize data. Assuming all sessions the same
#each must sum to 1
theta.marked.init <- c(0.5,0.25,0.25)
theta.unmarked.init <- c(0,0.75,0.25) #first cell for theta.unmarked.init should be 0

inits <- list(p0=rep(0.5,N.session),lam0=rep(1,N.session),sigma=rep(1,N.session)) #initializing with 1 parameter per session, just set all to same value

#This function structures the simulated data to fit the model in Nimble (some more restructing below)
nimbuild <- init.SMR.multisession.Dcov.Generalized(data,inits,M=M)
#plot to check s inits
for(g in 1:N.session){
  image(data[[g]]$x.vals,data[[g]]$y.vals,matrix(data[[g]]$D.cov*data[[g]]$InSS,data[[g]]$n.cells.x,data[[g]]$n.cells.y),
        main=paste("Session",g),xlab="X",ylab="Y",col=cols1)
  points(nimbuild$s[g,,],pch=16) #initialized activity centers
  points(X.sight[[g]],pch=4,col="lightblue")
  points(X.mark[[g]],pch=4)
  for(i in 1:n.marked[g]){
    trapcaps1 <- which(data[[g]]$y.mark[i,]>0)
    trapcaps2 <- which(data[[g]]$y.mID[i,]>0)
    traps <-  rbind(X.mark[[g]][trapcaps1,],X.sight[[g]][trapcaps2,])
    s <- nimbuild$s[g,i,]
    points(s[1],s[2],col="goldenrod",pch=16)
    if(nrow(traps)>0){
      for(j in 1:nrow(traps)){
        lines(x=c(s[1],traps[j,1]),y=c(s[2],traps[j,2]),col="goldenrod")
      }
    }
  }
}

#inits for nimble
N.init <- rowSums(nimbuild$z.init,na.rm=TRUE)
D0.init <- (rowSums(nimbuild$z.init))/(rowSums(nimbuild$InSS)*nimbuild$cellArea)

Niminits <- list(N=N.init,z=nimbuild$z.init,s=nimbuild$s.init,D0=D0.init,D.beta1=rep(0,N.session),
                 p0.fixed=0.5,lam0.fixed=1,sigma.fixed=1)

#constants for Nimble
# #If you do not have telemetry use these. Make sure telemetry BUGS code is commented out.
# constants <- list(N.session=N.session,n.marked=n.marked,M=M,J.mark=nimbuild$J.mark,J.sight=nimbuild$J.sight,
#                   K1D.mark=nimbuild$K1D.mark,K1D.sight=nimbuild$K1D.sight,
#                   xlim=nimbuild$xlim,ylim=nimbuild$ylim,D.cov=nimbuild$D.cov,res=nimbuild$res,
#                   cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells)
# 
# # Supply data to Nimble. marginalized data formatted in nimbuild object
# Nimdata <- list(y.mark=nimbuild$y.mark, #marking process
#                 y.mID=nimbuild$y.mID, #marked with ID
#                 y.mnoID=nimbuild$y.mnoID, #marked without ID
#                 y.um=nimbuild$y.um, #unmarked
#                 y.unk=nimbuild$y.unk, #unk marked status
#                 X.mark=nimbuild$X.mark,X.sight=nimbuild$X.sight,cells=nimbuild$cells,
#                 dummy.data=nimbuild$dummy.data,InSS=nimbuild$InSS)

#If you have telemetry use these instead. Make sure to uncomment telemetry BUGS code.
constants <- list(N.session=N.session,n.marked=n.marked,M=M,J.mark=nimbuild$J.mark,J.sight=nimbuild$J.sight,
                  K1D.mark=nimbuild$K1D.mark,K1D.sight=nimbuild$K1D.sight,
                  xlim=nimbuild$xlim,ylim=nimbuild$ylim,D.cov=nimbuild$D.cov,res=nimbuild$res,
                  cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells,
                  #telemetry stuff
                  tel.inds=nimbuild$tel.inds,
                  n.tel.inds=nimbuild$n.tel.inds,n.locs.ind=nimbuild$n.locs.ind)
Nimdata <- list(y.mark=nimbuild$y.mark, #marking process
                y.mID=nimbuild$y.mID, #marked with ID
                y.mnoID=nimbuild$y.mnoID, #marked without ID
                y.um=nimbuild$y.um, #unmarked
                y.unk=nimbuild$y.unk, #unk marked status
                X.mark=nimbuild$X.mark,X.sight=nimbuild$X.sight,cells=nimbuild$cells,
                dummy.data=nimbuild$dummy.data,InSS=nimbuild$InSS,locs=nimbuild$locs)

# set parameters to monitor
parameters <- c('D0','lambda.N','p0.fixed','lam0.fixed','sigma.fixed','theta.marked.fixed','theta.unmarked.fixed',
                'N','D.beta1')
parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting

# Build the model, configure the mcmc, and compile
start.time <- Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
#use block sampler below for 'D0.M','D0.UM','D.beta1'
config.nodes <- c('p0.fixed','lam0.fixed','sigma.fixed','theta.marked.fixed','theta.unmarked.fixed[2:3]')
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                      monitors2=parameters2,thin2=10,nodes=config.nodes)

# how many z proposals per iteration per session for marked and unmarked individuals?
z.ups <- round((M-n.marked)*0.25) #doing 25% of M-n.marked

for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  pd.nodes <- Rmodel$expandNodeNames(paste("pd[",g,",1:",M[g],",1:",J.mark[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",1:",M[g],",1:",J.sight[g],"]"))
  y.mark.nodes <- Rmodel$expandNodeNames(paste("y.mark[",g,",1:",M[g],",1:",J.mark[g],"]"))
  y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[",g,",1:",J.sight[g],"]"))
  y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[",g,",1:",J.sight[g],"]"))
  bigLam.unmarked.nodes <- Rmodel$expandNodeNames(paste("bigLam.unmarked[",g,",1:",J.sight[g],"]")) #only need this in calcNodes
  lam.um.nodes <- Rmodel$expandNodeNames(paste("lam.um[",g,",1:",J.sight[g],"]"))
  lam.unk.nodes <- Rmodel$expandNodeNames(paste("lam.unk[",g,",1:",J.sight[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",1:",M[g],"]"))
  calcNodes <- c(N.node,pd.nodes,lam.nodes,bigLam.unmarked.nodes,
                 lam.um.nodes,lam.unk.nodes,y.mark.nodes,y.um.nodes,y.unk.nodes)
  conf$addSampler(target = paste("N"),
                  type = 'zSampler',control = list(g=g,z.ups=z.ups[g],
                                                   J.sight=J.sight[g],J.mark=J.mark[g],
                                                   n.marked=n.marked[g],M=M[g],
                                                   pd.nodes=pd.nodes,lam.nodes=lam.nodes,
                                                   lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                   y.mark.nodes=y.mark.nodes,y.um.nodes=y.um.nodes,
                                                   y.unk.nodes=y.unk.nodes,
                                                   N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

#add sSampler
#if no telemetry,
# for(g in 1:N.session){
#   for(i in 1:M[g]){
#     conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
#                     type = 'sSampler',control=list(i=i,g=g,J.mark=J.mark[g],J.sight=J.sight[g],n.cells=nimbuild$n.cells[g],
#                                                    n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
#                                                    xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
#                                                    n.marked=n.marked[g],n.locs.ind=0,scale=1),silent = TRUE)
#     #scale parameter here is just the starting scale. It will be tuned.
#   }
# }
#if telemetry
for(g in 1:N.session){
  for(i in 1:M[g]){
    if(i %in% nimbuild$tel.inds[g,]){#inds with telemetry
      conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                      type = 'sSampler',control=list(i=i,g=g,J.mark=J.mark[g],J.sight=J.sight[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
                                                     n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                     xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
                                                     n.locs.ind=nimbuild$n.locs.ind[g,i],n.marked=n.marked[g],scale=1),silent = TRUE)
      #scale parameter here is just the starting scale. It will be tuned.
    }else{ #inds with no telemetry
      conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                      type = 'sSampler',control=list(i=i,g=g,J.mark=J.mark[g],J.sight=J.sight[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
                                                     n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                     xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
                                                     n.locs.ind=0,n.marked=n.marked[g],scale=1),silent = TRUE)
    }
  }
}

#can add block sampler if lam0, sigma, and/or lambda.N posteriors highly correlated
#probably not needed with telemetry
conf$addSampler(target = c("lam0.fixed","sigma.fixed"),type = 'RW_block',
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
unlist(lapply(data,function(x){x$lambda.N}))#expected Ns
unlist(lapply(data,function(x){x$N})) #realized Ns

tmp <- cor(mvSamples[200:nrow(mvSamples),])
diag(tmp) <- NA
which(abs(tmp)>0.5,arr.ind=TRUE)

#Important! If N[g] hits M[g] during sampling, raise M[g]. 

#Important! If N hits M during sampling, raise M. 

#Look at cell-level expected density estimates, compare to truth
mvSamples2  <-  as.matrix(Cmcmc$mvSamples2)
lambda.cell.idx <- grep("lambda.cell",colnames(mvSamples2))
n.cells.max <- max(n.cells)
lambda.cell.idx <- matrix(lambda.cell.idx,N.session,n.cells.max)
D0.idx <- grep("D0",colnames(mvSamples2))
burnin2 <- 10 #consider nt2 thinning rate when setting burnin2

#compare expected D plot to truth
#image will show posterior means
n.iter.use <- burnin2:nrow(mvSamples2)
lambda.cell.post <- array(NA,dim=c(N.session,n.cells.max,length(n.iter.use)))
lambda.cell <- lambda.cell.ests <- array(NA,dim=c(N.session,n.cells.max))
lambda.cell.HPDs <- array(NA,dim=c(N.session,n.cells.max,2))
for(g in 1:N.session){
  lambda.cell[g,1:n.cells[g]] <- exp(D.beta0[g] + D.beta1[g]*D.cov[[g]])*cellArea[g]
  lambda.cell.post[g,1:n.cells[g],] <- t(cellArea[g]*mvSamples2[n.iter.use,D0.idx[g]]*
                                           mvSamples2[n.iter.use,lambda.cell.idx[g,1:n.cells[g]]])
  lambda.cell.ests[g,1:n.cells[g]] <- rowMeans(lambda.cell.post[g,1:n.cells[g],])
  lambda.cell.HPDs[g,1:n.cells[g],] <- HPDinterval(mcmc(t(lambda.cell.post[g,1:n.cells[g],])))
  #remove nonhabitat (or not, comment out)
  lambda.cell[g,InSS[[g]]==0] <- NA
  lambda.cell.ests[g,InSS[[g]]==0] <- NA
}

par(mfrow=c(1,1),ask=FALSE)
for(g in 1:N.session){
  zlim <- range(c(lambda.cell[g,],lambda.cell.ests[g,]),na.rm=TRUE) #use same zlim for plots below
  #truth
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"True Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
  #estimate, posterior means
  image(x.vals[[g]],y.vals[[g]],matrix(lambda.cell.ests[g,1:n.cells[g]],n.cells.x[g],n.cells.y[g]),
        main=paste("Session",g,"Est Expected Density"),zlim=zlim)
  points(X[[g]],pch=4)
}
#cell ests and 95% HPDs vs. truth. 
#Need a lot of posterior samples for accurate 95% HPDs, if not, will look "jagged"
for(g in 1:N.session){
  idx <- order(lambda.cell[g,1:n.cells[g]])
  plot(lambda.cell.ests[g,1:n.cells[g]][idx]~lambda.cell[g,1:n.cells[g]][idx],type="l",lwd=2,
       main=paste("Session",g,"True vs. Estimated Density"),ylim=range(lambda.cell.HPDs[g,1:n.cells[g],]))
  lines(lambda.cell.HPDs[g,1:n.cells[g],1][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  lines(lambda.cell.HPDs[g,1:n.cells[g],2][idx]~lambda.cell[g,1:n.cells[g]][idx],lty=2)
  abline(0,1,col="darkred",lwd=2) #1:1 expectation
}
