args <- commandArgs(TRUE)
d <- as.numeric(args[1])
seed <- as.numeric(args[2])

setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Dcov_S1_fits")
filename <- paste0("S1_",d,".chain_1.RData")

files.done <- list.files()
if(!filename%in%files.done){#skip if already run
  
  nt <- 5
  n.iters <- 15000
  
  set.seed(seed)
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Dcov_S1_datasets")
  load(paste0("S1_data_",d,".RData"))
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/")
  library(nimble)
  source("sim.SMR.Dcov.R")
  source("NimbleModel SMR Poisson Dcov DA2 Marginal.R")
  source("NimbleFunctions SMR Poisson Dcov DA2 Marginal.R")
  source("init.SMR.Dcov.R")
  source("sSampler Poisson Dcov Marginal.R")
  source("mask.check.R")
  
  #If using Nimble version 0.13.1 and you must run this line
  nimbleOptions(determinePredictiveNodesInModel = FALSE)
  M <- 225
  X <- data$X
  J <- nrow(X)
  K <- data$K
  n.marked <- data$n.marked
  #Need some inits to initialize data
  #Use reasonable inits for lam0 and sigma since we check to make sure initial observation
  #model likelihood is finite
  #also use this function checks to make sure theta.marked and theta.unmarked inits are in
  #the correct structure. 
  inits <- list(lam0=1,sigma=1)
  
  #This function structures the simulated data to fit the model in Nimble (some more restructing below)
  #Also checks some inits
  nimbuild <- init.SMR.Dcov(data,inits,M=M)
  
  #inits for nimble
  D0.init <- (sum(nimbuild$z))/(sum(data$InSS)*data$res^2)
  
  #must initialize N.M and N.UM to be consistent with z. speeds converge to set consistent with lambda.N.M/UM
  Niminits <- list(z=nimbuild$z,s=nimbuild$s,D0=D0.init,D.beta1=0,
                   N=sum(nimbuild$z),
                   theta.unmarked=c(0,0.5,0.5),lam0=inits$lam0,sigma=inits$sigma)
  
  dummy.data <- rep(0,M) #dummy data not used, doesn't really matter what the values are
  
  #Use this if you do not have telemetry. Make sure telemetry commented out in model file
  constants <- list(n.marked=n.marked,M=M,J=J,K1D=data$K1D,
                    D.cov=data$D.cov,cellArea=data$cellArea,n.cells=data$n.cells,
                    xlim=nimbuild$xlim,ylim=nimbuild$ylim,res=data$res)
  
  Nimdata <- list(y.mID=nimbuild$y.mID, #marked with ID
                  y.mnoID=nimbuild$y.mnoID, #marked without ID
                  y.um=nimbuild$y.um, #unmarked
                  y.unk=nimbuild$y.unk, #unk marked status
                  X=as.matrix(X),
                  dummy.data=dummy.data,cells=data$cells,InSS=data$InSS)
  
  # set parameters to monitor
  parameters <- c('D0','lambda.N','lam0','sigma','theta.marked','theta.unmarked',
                  'N','D.beta1')
  parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting
  
  # Build the model, configure the mcmc, and compile
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
  #use block sampler below for 'D0.M','D0.UM','D.beta1'
  config.nodes <- c('lam0','sigma','theta.marked','theta.unmarked[2:3]')
  conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,
                        monitors2=parameters2,thin2=10,nodes=config.nodes)
  
  # how many z proposals per iteration per session for marked (if updated), unmarked?
  z.ups <- round(M*0.25) #doing 25% of M
  #nodes used for update, calcNodes + z nodes
  lam.nodes <- Rmodel$expandNodeNames("lam")
  y.um.nodes <- Rmodel$expandNodeNames("y.um")
  y.unk.nodes <- Rmodel$expandNodeNames("y.unk")
  bigLam.unmarked.nodes <- Rmodel$expandNodeNames("bigLam.unmarked") #only need this in calcNodes
  lam.um.nodes <- Rmodel$expandNodeNames("lam.um")
  lam.unk.nodes <- Rmodel$expandNodeNames("lam.unk")
  N.node <- Rmodel$expandNodeNames("N")
  z.nodes <- Rmodel$expandNodeNames("z")
  calcNodes <- c(N.node,lam.nodes,bigLam.unmarked.nodes,
                 lam.um.nodes,lam.unk.nodes,y.um.nodes,y.unk.nodes)
  conf$addSampler(target = paste("N"),
                  type = 'zSampler',control = list(z.ups=z.ups,J=J,K=K,
                                                   n.marked=n.marked,M=M,N.node=N.node,
                                                   y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                                   lam.nodes=lam.nodes,lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                   z.nodes=z.nodes,calcNodes=calcNodes),
                  silent = TRUE)
  
  #add sSampler
  # if no telemetry,
  for(i in 1:M){
    conf$addSampler(target = paste("s[",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,J=J,res=data$res,n.cells.x=data$n.cells.x,n.cells.y=data$n.cells.y,
                                                   xlim=nimbuild$xlim,ylim=nimbuild$ylim,
                                                   n.marked=n.marked,n.locs.ind=0,
                                                   scale=0.25),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }

  
  #can add block sampler if lam0, sigma, lambda.N.UM, and/or lambda.N.M (if N.M unknown) posteriors highly correlated
  #RW_block faster, AFslice slower, but mixes better
  conf$addSampler(target = c("lam0","sigma"),type = 'RW_block',
                  control = list(adaptive=TRUE),silent = TRUE)
  conf$addSampler(target = c("D0","D.beta1"),
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
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Dcov_S1_fits")
  save(out,file=filename)
}