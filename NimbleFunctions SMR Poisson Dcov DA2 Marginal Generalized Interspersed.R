dCell <- nimbleFunction(
  run = function(x = double(0), pi.cell = double(0),log = integer(0)) {
    returnType(double(0))
    logProb <- log(pi.cell)
    return(logProb)
  }
)

#make dummy random number generator to make nimble happy
rCell <- nimbleFunction(
  run = function(n = integer(0),pi.cell = double(0)) {
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
# Function to calculate detection rate, but skip when z=0
GetDetectionProb <- nimbleFunction(
  run = function(s = double(1), p0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- p0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)
#Vectorized observation model that also prevents z from being turned off if an unmarked ind currently has samples.
#also skips likelihood eval when z=0
dBinomialVector <- nimbleFunction(
  run = function(x = double(1), pd = double(1), K1D = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      if(sum(x)>0){ #need this so z is not turned off if samples allocated to individual
        return(-Inf)
      }else{
        return(0)
      }
    }else{
      logProb <- sum(dbinom(x, size = K1D, prob = pd, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rBinomialVector <- nimbleFunction(
  run = function(n = integer(0),pd = double(1), K1D = double(1), z = double(0)) {
    returnType(double(1))
    J <- nimDim(pd)[1]
    out <- numeric(J,value=0)
    return(out)
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
  setup = function(model, mvSaved, target, control){
    J.sight <- control$J.sight
    K.sight <- control$K.sight
    n.marked <- control$n.marked
    M <- control$M
    z.ups <- control$z.ups
    y.mark.nodes <- control$y.mark.nodes
    y.um.nodes <- control$y.um.nodes
    y.unk.nodes <- control$y.unk.nodes
    pd.nodes <- control$pd.nodes
    lam.nodes <- control$lam.nodes
    lam.um.nodes <- control$lam.um.nodes
    lam.unk.nodes <- control$lam.unk.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function(){
    bigLam.unmarked.initial <- model$bigLam.unmarked
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      if(updown==0){#subtract
        reject <- FALSE #we auto reject if you select a detected individual

        #find all z's currently on *including marked individuals*
        z.on <- which(model$z==1)
        n.z.on <- length(z.on)
        if(n.z.on>0){ #skip if no unmarked z's to turn off, otherwise nimble will crash
          pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
          pick <- z.on[pick]
          
          #prereject turning off any marked individuals or if there is a single unmarked individual
          if(model$N[1]==(n.marked+1)|pick<=n.marked){
            reject <- TRUE
          }
          if(!reject){
            #get initial logprobs for N and y
            lp.initial.N <- model$getLogProb(N.node)
            lp.initial.y.mark <- model$getLogProb(y.mark.nodes[pick])
            lp.initial.y.um <- model$getLogProb(y.um.nodes)
            lp.initial.y.unk <- model$getLogProb(y.unk.nodes)

            #propose new N/z
            model$N[1] <<-  model$N[1] - 1
            model$z[pick] <<- 0

            #turn off
            model$calculate(pd.nodes[pick])
            bigLam.unmarked.proposed <- bigLam.unmarked.initial - model$lam[pick,] #subtract these out before calculate
            model$calculate(lam.nodes[pick])
            model$bigLam.unmarked <<- bigLam.unmarked.proposed
            model$calculate(lam.um.nodes)
            model$calculate(lam.unk.nodes)

            #get proposed logprobs for N and y
            lp.proposed.N <- model$calculate(N.node)
            lp.proposed.y.mark <- model$calculate(y.mark.nodes[pick])
            lp.proposed.y.um <- model$calculate(y.um.nodes)
            lp.proposed.y.unk <- model$calculate(y.unk.nodes)

            #MH step
            log_MH_ratio <- (lp.proposed.N + lp.proposed.y.mark + lp.proposed.y.um + lp.proposed.y.unk) -
              (lp.initial.N + lp.initial.y.mark + lp.initial.y.um + lp.initial.y.unk)
            accept <- decide(log_MH_ratio)
            if(accept) {
              mvSaved["N",1][1] <<- model[["N"]]
              mvSaved["z",1][pick] <<- model[["z"]][pick]
              mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
              mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
              mvSaved["bigLam.unmarked",1][1:J.sight] <<- model[["bigLam.unmarked"]][1:J.sight]
              mvSaved["lam.um",1][1:J.sight,1:K.sight] <<- model[["lam.um"]][1:J.sight,1:K.sight]
              mvSaved["lam.unk",1][1:J.sight,1:K.sight] <<- model[["lam.unk"]][1:J.sight,1:K.sight]
              bigLam.unmarked.initial <- bigLam.unmarked.proposed
            }else{
              model[["N"]] <<- mvSaved["N",1][1]
              model[["z"]][pick] <<- mvSaved["z",1][pick]
              model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
              model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
              model[["bigLam.unmarked"]][1:J.sight] <<- mvSaved["bigLam.unmarked",1][1:J.sight]
              model[["lam.um"]][1:J.sight,1:K.sight] <<- mvSaved["lam.um",1][1:J.sight,1:K.sight]
              model[["lam.unk"]][1:J.sight,1:K.sight] <<- mvSaved["lam.unk",1][1:J.sight,1:K.sight]
              model$calculate(y.mark.nodes[pick])
              model$calculate(y.um.nodes)
              model$calculate(y.unk.nodes)
              model$calculate(N.node)
            }
          }
        }
      }else{#add
        if(model$N[1] < M){ #cannot update if z maxed out. Need to raise M

          #find all z's currently off.
          z.off <- which(model$z==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]
          
          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y.mark <- model$getLogProb(y.mark.nodes[pick])
          lp.initial.y.um <- model$getLogProb(y.um.nodes)
          lp.initial.y.unk <- model$getLogProb(y.unk.nodes)
          
          #propose new N/z
          model$N[1] <<-  model$N[1] + 1
          model$z[pick] <<- 1

          #turn on
          model$calculate(pd.nodes[pick])
          model$calculate(lam.nodes[pick])
          bigLam.unmarked.proposed <- bigLam.unmarked.initial + model$lam[pick,] #add these in after calculate
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          model$calculate(lam.um.nodes)
          model$calculate(lam.unk.nodes)

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y.mark <- model$calculate(y.mark.nodes[pick])
          lp.proposed.y.um <- model$calculate(y.um.nodes)
          lp.proposed.y.unk <- model$calculate(y.unk.nodes)

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y.mark + lp.proposed.y.um + lp.proposed.y.unk) -
            (lp.initial.N + lp.initial.y.mark + lp.initial.y.um + lp.initial.y.unk)
          accept <- decide(log_MH_ratio)
          if(accept){
            mvSaved["N",1][1] <<- model[["N"]]
            mvSaved["z",1][pick] <<- model[["z"]][pick]
            mvSaved["pd",1][pick,] <<- model[["pd"]][pick,]
            mvSaved["lam",1][pick,] <<- model[["lam"]][pick,]
            mvSaved["bigLam.unmarked",1][1:J.sight] <<- model[["bigLam.unmarked"]][1:J.sight]
            mvSaved["lam.um",1][1:J.sight,1:K.sight] <<- model[["lam.um"]][1:J.sight,1:K.sight]
            mvSaved["lam.unk",1][1:J.sight,1:K.sight] <<- model[["lam.unk"]][1:J.sight,1:K.sight]
            bigLam.unmarked.initial <- bigLam.unmarked.proposed
          }else{
            model[["N"]] <<- mvSaved["N",1][1]
            model[["z"]][pick] <<- mvSaved["z",1][pick]
            model[["pd"]][pick,] <<- mvSaved["pd",1][pick,]
            model[["lam"]][pick,] <<- mvSaved["lam",1][pick,]
            model[["bigLam.unmarked"]][1:J.sight] <<- mvSaved["bigLam.unmarked",1][1:J.sight]
            model[["lam.um"]][1:J.sight,1:K.sight] <<- mvSaved["lam.um",1][1:J.sight,1:K.sight]
            model[["lam.unk"]][1:J.sight,1:K.sight] <<- mvSaved["lam.unk",1][1:J.sight,1:K.sight]
            model$calculate(y.mark.nodes[pick])
            model$calculate(y.um.nodes)
            model$calculate(y.unk.nodes)
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