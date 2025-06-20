args <- commandArgs(TRUE)
d <- as.numeric(args[1])
seed <- as.numeric(args[2])

setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Multisession_tel_S1_fits")
filename <- paste0("S1_",d,".chain_1.RData")

files.done <- list.files()
if(!filename%in%files.done){#skip if already run
  
  nt <- 5
  n.iters <- 15000
  set.seed(seed)
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Multisession_tel_S1_datasets")
  load(paste0("S1_data_",d,".RData"))
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/")
  library(nimble)
  source("sim.SMR.multisession.Dcov.R")
  source("sim.SMR.Dcov.R") #used by multisession simulator
  source("NimbleModel SMR Multisession Poisson Dcov Marginal.R")
  source("NimbleFunctions SMR Multisession Poisson Dcov Marginal.R")
  source("init.SMR.multisession.Dcov.R")
  source("init.SMR.Dcov.R") #used by multisession initializer
  source("sSampler Multisession Poisson Dcov Marginal.R")
  source("mask.check.R")
  
  #If using Nimble version 0.13.1 and you must run this line
  nimbleOptions(determinePredictiveNodesInModel = FALSE)
  
  ####Fit model in Nimble####
  N.session <- 3
  M <- c(225,225,225) #Augmentation level for each session
  X <- lapply(data,function(x){x$X})
  J <- sapply(X,nrow)
  n.marked <- sapply(data,function(x){x$n.marked})
  
  #Need some inits to initialize data
  #Use reasonable inits for lam0 and sigma since we check to make sure initial observation
  #model likelihood is finite
  #provide inits for theta.marked and theta.unmarked to initialize data. Assuming all sessions the same
  #each must sum to 1
  theta.marked.init <- c(0.5,0.25,0.25)
  theta.unmarked.init <- c(0,0.75,0.25) #first cell for theta.unmarked.init should be 0
  
  inits <- list(lam0=rep(1,N.session),sigma=rep(1,N.session)) #initializing with 1 parameter per session, just set all to same value
  
  #This function structures the simulated data to fit the model in Nimble (some more restructing below)
  nimbuild <- init.SMR.multisession.Dcov(data,inits,M=M)
  
  #inits for nimble
  N.init <- rowSums(nimbuild$z.init,na.rm=TRUE)
  D0.init <- (rowSums(nimbuild$z.init))/(rowSums(nimbuild$InSS)*nimbuild$cellArea)
  
  Niminits <- list(N=N.init,z=nimbuild$z.init,s=nimbuild$s.init,D0=D0.init,D.beta1=rep(0,N.session),
                   lam0.fixed=1,sigma.fixed=1)
  
  constants <- list(N.session=N.session,n.marked=n.marked,M=M,J=J,K1D=nimbuild$K1D,
                    xlim=nimbuild$xlim,ylim=nimbuild$ylim,D.cov=nimbuild$D.cov,res=nimbuild$res,
                    cellArea=nimbuild$cellArea,n.cells=nimbuild$n.cells,
                    #telemetry stuff
                    tel.inds=nimbuild$tel.inds,
                    n.tel.inds=nimbuild$n.tel.inds,n.locs.ind=nimbuild$n.locs.ind)
  Nimdata <- list(y.mID=nimbuild$y.mID, #marked with ID
                  y.mnoID=nimbuild$y.mnoID, #marked without ID
                  y.um=nimbuild$y.um, #unmarked
                  y.unk=nimbuild$y.unk, #unk marked status
                  X=nimbuild$X,cells=nimbuild$cells,
                  dummy.data=nimbuild$dummy.data,InSS=nimbuild$InSS,locs=nimbuild$locs)
  
  # set parameters to monitor
  parameters <- c('D0','D.beta1','lambda.N','N','lam0.fixed','sigma.fixed',
                  'theta.marked.fixed','theta.unmarked.fixed')
  parameters2 <- c("lambda.cell","s.cell",'D0') #record D0 here for plotting
  nt <- 1 #thinning rate for parameters
  nt2 <- 10 #thinning rate for parameters2
  
  
  # Build the model, configure the mcmc, and compile
  start.time <- Sys.time()
  Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,inits=Niminits)
  #if you change the model file, make sure you make necessary changes here
  config.nodes <- c('lam0.fixed','sigma.fixed','theta.marked.fixed','theta.unmarked.fixed[2:3]')
  #this one for session-specific theta.marked and theta.unmarked
  # config.nodes <- c('lam0.fixed','sigma.fixed','theta.marked',paste('theta.unmarked[1:',N.session,',2:3',']'))
  conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,monitors2=parameters2, thin2=nt2,nodes=config.nodes)
  
  # how many z proposals per iteration per session for marked and unmarked individuals?
  z.ups <- round((M-n.marked)*0.25) #doing 25% of M-n.marked
  
  for(g in 1:N.session){
    #nodes used for update, calcNodes + z nodes
    lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",1:",M[g],",1:",J[g],"]"))
    y.um.nodes <- Rmodel$expandNodeNames(paste("y.um[",g,",1:",J[g],"]"))
    y.unk.nodes <- Rmodel$expandNodeNames(paste("y.unk[",g,",1:",J[g],"]"))
    bigLam.unmarked.nodes <- Rmodel$expandNodeNames(paste("bigLam.unmarked[",g,",1:",J[g],"]")) #only need this in calcNodes
    lam.um.nodes <- Rmodel$expandNodeNames(paste("lam.um[",g,",1:",J[g],"]"))
    lam.unk.nodes <- Rmodel$expandNodeNames(paste("lam.unk[",g,",1:",J[g],"]"))
    N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
    z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",1:",M[g],"]"))
    calcNodes <- c(N.node=N.node,lam.nodes,bigLam.unmarked.nodes,
                   lam.um.nodes,lam.unk.nodes,y.um.nodes,y.unk.nodes)
    conf$addSampler(target = paste("N"),
                    type = 'zSampler',control = list(g=g,z.ups=z.ups[g],
                                                     J=J[g],n.marked=n.marked[g],M=M[g],
                                                     y.um.nodes=y.um.nodes,y.unk.nodes=y.unk.nodes,
                                                     lam.nodes=lam.nodes,
                                                     lam.um.nodes=lam.um.nodes,lam.unk.nodes=lam.unk.nodes,
                                                     N.node=N.node,z.nodes=z.nodes,
                                                     calcNodes=calcNodes),
                    silent = TRUE)
  }
  
  #add sSampler
  for(g in 1:N.session){
    for(i in 1:M[g]){
      if(i %in% nimbuild$tel.inds[g,]){#inds with telemetry
        conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                        type = 'sSampler',control=list(i=i,g=g,J=J[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
                                                       n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                       xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
                                                       n.locs.ind=nimbuild$n.locs.ind[g,i],n.marked=n.marked[g],scale=1),silent = TRUE)
        #scale parameter here is just the starting scale. It will be tuned.
      }else{ #inds with no telemetry
        conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                        type = 'sSampler',control=list(i=i,g=g,J=J[g],res=nimbuild$res[g],n.cells=nimbuild$n.cells[g],
                                                       n.cells.x=nimbuild$n.cells.x[g],n.cells.y=nimbuild$n.cells.y[g],
                                                       xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],res=nimbuild$res[g],
                                                       n.locs.ind=0,n.marked=n.marked[g],scale=1),silent = TRUE)
      }
    }
  }
  
  #mixing seems generally good if you keep independent samplers and add block sampler (if data sparse enough that there is posterior correlation)
  # conf$removeSampler(c("lam0.fixed","sigma.fixed"))
  conf$addSampler(target = c("lam0.fixed","sigma.fixed"),type = 'RW_block',
                  control = list(adaptive=TRUE),silent = TRUE)
  #AF_slice pretty efficient for Dcov parameters. Block by session
  for(g in 1:N.session){
    conf$addSampler(target = c(paste("D0[",g,"]"),paste("D.beta1[",g,"]")),
                    type = 'AF_slice',control=list(adaptive=TRUE),silent = TRUE)
  }
  
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
  
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Multisession_tel_S1_fits")
  save(out,file=filename)
}