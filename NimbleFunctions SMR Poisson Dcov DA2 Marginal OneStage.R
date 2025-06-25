dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0), InSS = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(InSS==1){
      logProb <- log(pi.cell)
    }else{
      logProb <- -Inf
    }
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0), InSS = double(0)) {
    returnType(double(0))
    return(0)
  }
)

# Function to calculate detection rate, but skip when z=0
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)
#Vectorized observation model that also prevents z from being turned off if an unmarked ind currently has samples.
#also skips likelihood eval when z=0
dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lambda = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dpois(x, lambda = lambda, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lambda = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(lambda)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)
GetbigLam <- nimbleFunction(
  run = function(lam = double(2), z = double(1)){ 
    returnType(double(1))
    M <- nimDim(lam)[1]
    J <- nimDim(lam)[2]
    bigLam <- rep(0,J)
    for(i in 1:M){
      if(z[i]==1){
        bigLam <- bigLam + lam[i,]
      }
    }
    return(bigLam)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    J <- control$J
    n.marked <- control$n.marked
    M <- control$M
    z.ups <- control$z.ups
    y.all.nodes <- control$y.all.nodes
    lam.nodes <- control$lam.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    bigLam.all.initial <- model$bigLam.all
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        reject <- FALSE #we auto reject if you select a detected individual
        #find all z's currently on *excluding marked individuals*
        z.on <- which(model$z[(n.marked+1):M]==1) + n.marked
        n.z.on <- length(z.on)
        if(n.z.on>0){ #skip if no unmarked z's to turn off, otherwise nimble will crash
          pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
          pick <- z.on[pick]
          
          #prereject turning off all unmarked individuals
          if(model$N[1]==(n.marked+1)){ #is this the last unmarked individual?
            reject <- TRUE
          }
          if(!reject){
            #get initial logprobs for N and y
            lp.initial.N <- model$getLogProb(N.node)
            lp.initial.y.all <- model$getLogProb(y.all.nodes)
            
            #propose new N/z
            model$N[1] <<-  model$N[1] - 1
            model$z[pick] <<- 0
            
            #turn off
            bigLam.all.proposed <- bigLam.all.initial - model$lam[pick,] #subtract these out before calculate
            model$calculate(lam.nodes[pick])
            model$bigLam.all <<- bigLam.all.proposed
            
            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.node)
            lp.proposed.y.all <- model$calculate(y.all.nodes)
            
            #MH step
            log_MH_ratio <- (lp.proposed.N + lp.proposed.y.all) - (lp.initial.N + lp.initial.y.all)
            accept <- decide(log_MH_ratio)
            if(accept) {
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["z",1][pick] <<- model[["z"]][pick]
              mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
              mvSaved["bigLam.all",1][1:J] <<- model[["bigLam.all"]][1:J]
              bigLam.all.initial <- bigLam.all.proposed
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["z"]][pick] <<- mvSaved["z",1][pick]
              model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
              model[["bigLam.all"]][1:J] <<- mvSaved["bigLam.all",1][1:J]
              model$calculate(y.all.nodes)
              model$calculate(N.node)
            }
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M
          
          #find all z's currently off
          z.off <- which(model$z[1:M]==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y.all <- model$getLogProb(y.all.nodes)
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1
          
          #turn on
          model$calculate(lam.nodes[pick])
          bigLam.all.proposed <- bigLam.all.initial + model$lam[pick,] #add these in after calculate
          model$bigLam.all <<- bigLam.all.proposed
          
          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y.all <- model$calculate(y.all.nodes)
          
          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y.all) - (lp.initial.N + lp.initial.y.all)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["bigLam.all",1][1:J] <<- model[["bigLam.all"]][1:J]
            bigLam.all.initial <- bigLam.all.proposed
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["bigLam.all"]][1:J] <<- mvSaved["bigLam.all",1][1:J]
            model$calculate(y.all.nodes)
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)