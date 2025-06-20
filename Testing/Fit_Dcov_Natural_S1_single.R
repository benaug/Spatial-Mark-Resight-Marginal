args <- commandArgs(TRUE)
d <- as.numeric(args[1])
seed <- as.numeric(args[2])

setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Dcov_Natural_S1_fits")
filename <- paste0("S1_",d,".chain_1.RData")

files.done <- list.files()
if(!filename%in%files.done){#skip if already run
  
  nt <- 5
  n.iters <- 15000
  
  set.seed(seed)
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Dcov_Natural_S1_datasets")
  load(paste0("S1_data_",d,".RData"))
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/")
  library(nimble)
  source("sim.SMR.Dcov.Natural.R")
  source("NimbleModel SMR Poisson Dcov DA2 Marginal Natural.R")
  source("NimbleFunctions SMR Poisson Dcov DA2 Marginal Natural.R")
  source("init.SMR.Dcov.Natural.R")
  source("sSampler Poisson Dcov Marginal Natural.R")
  source("mask.check.R")
  
  #If using Nimble version 0.13.1 and you must run this line
  nimbleOptions(determinePredictiveNodesInModel = FALSE)
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
  # image(data$x.vals,data$y.vals,matrix(data$D.cov,data$n.cells.x,data$n.cells.y),main="D.cov",xlab="X",ylab="Y",col=cols1)
  # points(X,pch=4)
  # points(nimbuild$s,pch=16) #initialized activity centers
  
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
  Cmcmc$run(n.iters,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
  end.time <- Sys.time()
  end.time-start.time  # total time for compilation, replacing samplers, and fitting
  end.time-start.time2 # post-compilation run time
  
  mvSamples = as.matrix(Cmcmc$mvSamples)
  
  # burnin <- 50
  # plot(coda::mcmc(mvSamples[-c(1:burnin),]))
  
  out <- list(mvSamples=mvSamples,seed=seed,time=end.time-start.time2)
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Dcov_Natural_S1_fits")
  save(out,file=filename)
}