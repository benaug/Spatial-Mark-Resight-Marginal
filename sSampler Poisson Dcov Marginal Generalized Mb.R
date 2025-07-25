sSampler <- nimbleFunction(
  # name = 'sampler_RW',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    i <- control$i
    J.mark <- control$J.mark
    J.sight <- control$J.sight
    res <- control$res
    xlim <- control$xlim
    ylim <- control$ylim
    n.cells.x <- control$n.cells.x
    n.cells.y <- control$n.cells.y
    n.locs.ind <- control$n.locs.ind
    n.marked <- control$n.marked
    ## control list extraction
    # logScale            <- extractControlElement(control, 'log',                 FALSE)
    # reflective          <- extractControlElement(control, 'reflective',          FALSE)
    adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
    adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
    adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
    scale               <- extractControlElement(control, 'scale',               1)
    ## node list generation
    # targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    calcNodes <- model$getDependencies(target)
    s.nodes <- c(model$expandNodeNames(paste("s[",i,",",1:2,"]")),
                 model$expandNodeNames(paste("s.cell[",i,"]")),
                 model$expandNodeNames(paste("dummy.data[",i,"]")))
    #if we have telemetry for this individual, add locs to s.nodes
    if(n.locs.ind>0){
      loc.nodes <- model$expandNodeNames(paste("locs[",i,",1:",n.locs.ind,",",1:2,"]"))
      s.nodes <- c(s.nodes,loc.nodes)
    }
    y.mark.p.nodes <- model$expandNodeNames(paste("y.mark.p[",i,",1:",J.mark,"]"))
    y.mark.c.nodes <- model$expandNodeNames(paste("y.mark.c[",i,",1:",J.mark,"]"))
    y.mID.nodes <- model$expandNodeNames(paste("y.mID[",i,",1:",J.sight,"]"))
    y.mnoID.nodes <- model$expandNodeNames("y.mnoID")
    y.um.nodes <- model$expandNodeNames("y.um")
    y.unk.nodes <- model$expandNodeNames("y.unk")
    kern.nodes <- model$expandNodeNames(paste("kern[",i,",1:",J.mark,"]"))
    pd.p.nodes <- model$expandNodeNames(paste("pd.p[",i,",1:",J.mark,"]"))
    pd.c.nodes <- model$expandNodeNames(paste("pd.c[",i,",1:",J.mark,"]"))
    lam.nodes <- model$expandNodeNames(paste("lam[",i,",1:",J.sight,"]"))
    lam.mnoID.nodes <- model$expandNodeNames(paste("lam.mnoID[1:",J.sight,"]"))
    lam.um.nodes <- model$expandNodeNames(paste("lam.um[1:",J.sight,"]"))
    lam.unk.nodes <- model$expandNodeNames(paste("lam.unk[1:",J.sight,"]"))
    # calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
    # isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
    # calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
    # calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
    ## numeric value generation
    scaleOriginal <- scale
    timesRan      <- 0
    timesAccepted <- 0
    timesAdapted  <- 0
    scaleHistory  <- c(0, 0)   ## scaleHistory
    acceptanceHistory  <- c(0, 0)   ## scaleHistory
    if(nimbleOptions('MCMCsaveHistory')) {
      saveMCMChistory <- TRUE
    } else saveMCMChistory <- FALSE
    optimalAR     <- 0.44
    gamma1        <- 0
    ## checks
    # if(length(targetAsScalar) > 1)   stop('cannot use RW sampler on more than one target; try RW_block sampler')
    # if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
    # if(logScale & reflective)        stop('cannot use reflective RW sampler on a log scale (i.e. with options log=TRUE and reflective=TRUE')
    if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
    if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
  },
  run = function() {
    z <- model$z[i]
    if(z==0){ #propose from prior
      # #propose new cell
      model$s.cell[i] <<- rcat(1,model$pi.cell)
      #propose x and y in new cell
      s.cell.x <- model$s.cell[i]%%n.cells.x
      s.cell.y <- floor(model$s.cell[i]/n.cells.x)+1
      if(s.cell.x==0){
        s.cell.x <- n.cells.x
        s.cell.y <- s.cell.y-1
      }
      xlim.cell <- c(s.cell.x-1,s.cell.x)*res
      ylim.cell <- c(s.cell.y-1,s.cell.y)*res
      model$s[i,1:2] <<- c(runif(1, xlim.cell[1], xlim.cell[2]), runif(1, ylim.cell[1], ylim.cell[2]))
      model$calculate(s.nodes)
      copy(from = model, to = mvSaved, row = 1, nodes = s.nodes, logProb = TRUE)
    }else{ #MH
      s.cand <- c(rnorm(1,model$s[i,1],scale), rnorm(1,model$s[i,2],scale))
      inbox <- s.cand[1]< xlim[2] & s.cand[1]> xlim[1] & s.cand[2] < ylim[2] & s.cand[2] > ylim[1]
      if(inbox){
        #get initial logprobs
        lp_initial_s <- model$getLogProb(s.nodes)
        lp_initial_y.mark.p <- model$getLogProb(y.mark.p.nodes)
        lp_initial_y.mark.c <- model$getLogProb(y.mark.c.nodes)
        lp_initial_y.unk <- model$getLogProb(y.unk.nodes)
        if(i<=n.marked){ #if marked
          lp_initial_y.mID <- model$getLogProb(y.mID.nodes)
          lp_initial_y.mnoID <- model$getLogProb(y.mnoID.nodes)
          #pull this out of model object
          bigLam.marked.initial <- model$bigLam.marked
          #update proposed s
          model$s[i, 1:2] <<- s.cand
          lp_proposed_s <- model$calculate(s.nodes) #proposed logprob for s.nodes
          model$calculate(kern.nodes) #update kern nodes
          model$calculate(pd.p.nodes) #update pd nodes
          model$calculate(pd.c.nodes) #update pd nodes
          #subtract these out before calculating lam
          bigLam.marked.proposed <- bigLam.marked.initial - model$lam[i,]
          model$calculate(lam.nodes) #update lam nodes
          #add these in after calculating lam
          bigLam.marked.proposed <- bigLam.marked.proposed + model$lam[i,]
          #put bigLam in model object
          model$bigLam.marked <<- bigLam.marked.proposed
          model$calculate(lam.mnoID.nodes) #update after bigLam
          model$calculate(lam.unk.nodes) #update after bigLam
          lp_proposed_y.mark.p <- model$calculate(y.mark.p.nodes)
          lp_proposed_y.mark.c <- model$calculate(y.mark.c.nodes)
          lp_proposed_y.mID <- model$calculate(y.mID.nodes)
          lp_proposed_y.mnoID <- model$calculate(y.mnoID.nodes)
          lp_proposed_y.unk <- model$calculate(y.unk.nodes)
          lp_initial <- lp_initial_s + lp_initial_y.mark.p + lp_initial_y.mark.c + lp_initial_y.mID + lp_initial_y.mnoID + lp_initial_y.unk
          lp_proposed <- lp_proposed_s + lp_proposed_y.mark.p + lp_proposed_y.mark.c + lp_proposed_y.mID + lp_proposed_y.mnoID + lp_proposed_y.unk
        }else{ #else unmarked
          lp_initial_y.um <- model$getLogProb(y.um.nodes)
          #pull this out of model object
          bigLam.unmarked.initial <- model$bigLam.unmarked
          #update proposed s
          model$s[i, 1:2] <<- s.cand
          lp_proposed_s <- model$calculate(s.nodes) #proposed logprob for s.nodes
          model$calculate(kern.nodes) #update kern nodes
          model$calculate(pd.p.nodes) #update pd nodes
          model$calculate(pd.c.nodes) #update pd nodes
          #subtract these out before calculating lam
          bigLam.unmarked.proposed <- bigLam.unmarked.initial - model$lam[i,]
          model$calculate(lam.nodes) #update lam nodes
          #add these in after calculating lam
          bigLam.unmarked.proposed <- bigLam.unmarked.proposed + model$lam[i,]
          #put bigLam in model object
          model$bigLam.unmarked <<- bigLam.unmarked.proposed
          model$calculate(lam.um.nodes) #update after bigLam
          model$calculate(lam.unk.nodes) #update after bigLam
          lp_proposed_y.mark.p <- model$calculate(y.mark.p.nodes)
          lp_proposed_y.mark.c <- model$calculate(y.mark.c.nodes)
          lp_proposed_y.um <- model$calculate(y.um.nodes)
          lp_proposed_y.unk <- model$calculate(y.unk.nodes)
          lp_initial <- lp_initial_s + lp_initial_y.mark.p + lp_initial_y.mark.c + lp_initial_y.um + lp_initial_y.unk
          lp_proposed <- lp_proposed_s + lp_proposed_y.mark.p + lp_proposed_y.mark.c + lp_proposed_y.um + lp_proposed_y.unk
        }
        log_MH_ratio <- lp_proposed - lp_initial
        accept <- decide(log_MH_ratio)
        if(accept) {
          copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
        } else {
          copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
        }
        if(adaptive){ #we only tune for z=0 proposals
          adaptiveProcedure(accept)
        }
      }
    }
  },
  methods = list(
    adaptiveProcedure = function(jump = logical()) {
      timesRan <<- timesRan + 1
      if(jump)     timesAccepted <<- timesAccepted + 1
      if(timesRan %% adaptInterval == 0) {
        acceptanceRate <- timesAccepted / timesRan
        timesAdapted <<- timesAdapted + 1
        if(saveMCMChistory) {
          setSize(scaleHistory, timesAdapted)                 ## scaleHistory
          scaleHistory[timesAdapted] <<- scale                ## scaleHistory
          setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
          acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
        }
        gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
        gamma2 <- 10 * gamma1
        adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
        scale <<- scale * adaptFactor
        ## If there are upper and lower bounds, enforce a maximum scale of
        ## 0.5 * (upper-lower).  This is arbitrary but reasonable.
        ## Otherwise, for a poorly-informed posterior,
        ## the scale could grow without bound to try to reduce
        ## acceptance probability.  This creates enormous cost of
        ## reflections.
        # if(reflective) {
        #   lower <- model$getBound(target, 'lower')
        #   upper <- model$getBound(target, 'upper')
        #   if(scale >= 0.5*(upper-lower)) {
        #     scale <<- 0.5*(upper-lower)
        #   }
        # }
        timesRan <<- 0
        timesAccepted <<- 0
      }
    },
    getScaleHistory = function() {       ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(scaleHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },          
    getAcceptanceHistory = function() {  ## scaleHistory
      returnType(double(1))
      if(saveMCMChistory) {
        return(acceptanceHistory)
      } else {
        print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
        return(numeric(1, 0))
      }
    },
    ##getScaleHistoryExpanded = function() {                                                 ## scaleHistory
    ##    scaleHistoryExpanded <- numeric(timesAdapted*adaptInterval, init=FALSE)            ## scaleHistory
    ##    for(iTA in 1:timesAdapted)                                                         ## scaleHistory
    ##        for(j in 1:adaptInterval)                                                      ## scaleHistory
    ##            scaleHistoryExpanded[(iTA-1)*adaptInterval+j] <- scaleHistory[iTA]         ## scaleHistory
    ##    returnType(double(1)); return(scaleHistoryExpanded) },                             ## scaleHistory
    reset = function() {
      scale <<- scaleOriginal
      timesRan      <<- 0
      timesAccepted <<- 0
      timesAdapted  <<- 0
      if(saveMCMChistory) {
        scaleHistory  <<- c(0, 0)    ## scaleHistory
        acceptanceHistory  <<- c(0, 0)
      }
      gamma1 <<- 0
    }
  )
)