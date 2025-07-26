NimModel <- nimbleCode({
  #Density covariates, marked and unmarked have separate intercepts
  #baseline density
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D.beta1 ~ dnorm(0,sd=10)
  
  #detection function priors
  lam0 ~ dunif(0,15)
  sigma ~ dunif(0,10)

  #sample type observation model priors (Dirichlet)
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])

  #Density model
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space.
  
  #Marked individuals first
  for(i in 1:n.marked){
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
    for(k in 1:K){
      #marked and identified detections. only count when marked.status[i,k] = 1
      y.mID[i,1:J,k] ~ dPoissonVector(lam[i,1:J]*K2D[1:J,k]*marked.status[i,k]*theta.marked[1],z=z[i])
    }
  }
  
  #Then unmarked individuals
  for(i in (n.marked+1):M){
    #dunif() here implies uniform distribution within a grid cell
    #also tells nimble s's are in continuous space, not discrete
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]],InSS=InSS[s.cell[i]])
    lam[i,1:J] <- GetDetectionRate(s = s[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=z[i])
  }#custom Metropolis-Hastings update for N.UM/z[(n.marked+1):M]
  
  #Unidentified detections by type
  #1 marked with no ID detections
  for(k in 1:K){
    bigLam.marked[1:J,k] <- GetbigLam(lam=lam[1:n.marked,1:J],z=z[1:n.marked]*marked.status[1:n.marked,k])
    lam.mnoID[1:J,k] <- bigLam.marked[1:J,k]*K2D[1:J,k]*theta.marked[2]
    y.mnoID[1:J,k] ~ dPoissonVector(lam.mnoID[1:J,k],z=1) #plug in z=1 to reuse dPoissonVector
  }
  
  #2 unmarked detections
  bigLam.unmarked[1:J] <- GetbigLam(lam=lam[(n.marked+1):M,1:J],z=z[(n.marked+1):M])
  for(k in 1:K){
    bigLam.premarked[1:J,k] <- GetbigLam(lam=lam[1:n.marked,1:J],z=z[1:n.marked]*(1-marked.status[1:n.marked,k]))
    lam.um[1:J,k] <- bigLam.unmarked[1:J]*K2D[1:J,k]*theta.unmarked[2] + #always unmarked individuals
      bigLam.premarked[1:J,k]*K2D[1:J,k]*theta.unmarked[2] #marked inds when not marked
    y.um[1:J,k] ~ dPoissonVector(lam.um[1:J,k],z=1) #plug in z=1 to reuse dPoissonVector
  }
  
  #3 unknown marked status
  for(k in 1:K){
    lam.unk[1:J,k] <- bigLam.marked[1:J,k]*K2D[1:J,k]*theta.marked[3] + #marked individuals when marked
      bigLam.premarked[1:J,k]*K2D[1:J,k]*theta.unmarked[3] + #marked individuals when not marked
      bigLam.unmarked[1:J]*K2D[1:J,k]*theta.unmarked[3] #always unmarked individuals
    y.unk[1:J,k] ~ dPoissonVector(lam.unk[1:J,k],z=1) #plug in z=1 to reuse dPoissonVector
  }
  
  #If you have telemetry
  for(i in 1:n.tel.inds){
    for(m in 1:n.locs.ind[i]){
      locs[tel.inds[i],m,1] ~ dnorm(s[tel.inds[i],1],sd=sigma)
      locs[tel.inds[i],m,2] ~ dnorm(s[tel.inds[i],2],sd=sigma)
    }
  }
})# end model