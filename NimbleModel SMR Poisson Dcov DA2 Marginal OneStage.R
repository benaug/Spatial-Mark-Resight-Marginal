NimModel <- nimbleCode({
  #baseline density
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  D.beta1 ~ dnorm(0,sd=10)
  # D.beta0 ~ dnorm(0,sd=10)

  #detection function priors
  lam0 ~ dunif(0,15)
  sigma ~ dunif(0,10)

  #marked individual ID thinning rate, equivalent to theta.marked[1]
  theta.thin ~ dunif(0,1)

  #Density model
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space.
  
  #1) marked individuals
  for(i in 1:n.marked){
    #uniform assumption here
    s.m[i,1] ~ dunif(xlim[1],xlim[2])
    s.m[i,2] ~ dunif(ylim[1],ylim[2])
    lam.m[i,1:J] <- GetDetectionRate(s = s.m[i,1:2], X = X[1:J,1:2], J=J,sigma=sigma, lam0=lam0, z=1)
    y.mID[i,1:J] ~ dPoissonVector(lam.m[i,1:J]*K1D[1:J]*theta.thin,z=z.m[i]) #marked and identified detections
  }
  #marked with no ID detections
  bigLam.marked[1:J] <- GetbigLam(lam=lam.m[1:n.marked,1:J],z=z.m[1:n.marked])
  lam.mnoID[1:J] <- bigLam.marked[1:J]*K1D[1:J]*(1-theta.thin)
  y.mnoID[1:J] ~ dPoissonVector(lam.mnoID[1:J],z=1) #plug in z=1 to reuse dPoissonVector
  
  #2) all individuals treated as unmarked
  for(i in 1:M){
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
  }#custom Metropolis-Hastings update for N/z
  bigLam.all[1:J] <- GetbigLam(lam=lam[1:M,1:J],z=z[1:M])
  y.all[1:J] ~ dPoissonVector(bigLam.all[1:J]*K1D[1:J],z=1) #plug in z=1 to reuse dPoissonVector
  
  #If you have telemetry
  for(i in 1:n.tel.inds){
    for(m in 1:n.locs.ind[i]){
      locs[tel.inds[i],m,1] ~ dnorm(s.m[tel.inds[i],1],sd=sigma)
      locs[tel.inds[i],m,2] ~ dnorm(s.m[tel.inds[i],2],sd=sigma)
    }
  }
})# end model