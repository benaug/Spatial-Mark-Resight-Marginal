NimModel <- nimbleCode({
  #detection function priors - shared across sessions
  lam0.fixed ~ dunif(0,15)
  sigma.fixed ~ dunif(0,10)
  
  #marked individual ID thinning rate, equivalent to theta.marked[1]
  theta.thin.fixed ~ dunif(0,1)
  
  for(g in 1:N.session){
    #not sharing Dcov parameters over sessions, but can do that
    D0[g] ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
    D.beta1[g] ~ dnorm(0,sd=10)
    #Density model
    D.intercept[g] <- D0[g]*cellArea[g]
    lambda.cell[g,1:n.cells[g]] <- InSS[g,1:n.cells[g]]*exp(D.beta1[g]*D.cov[g,1:n.cells[g]])
    pi.cell[g,1:n.cells[g]] <- lambda.cell[g,1:n.cells[g]]/pi.denom[g] #expected proportion of total N in cell c
    pi.denom[g] <- sum(lambda.cell[g,1:n.cells[g]])
    lambda.N[g] <- D.intercept[g]*pi.denom[g] #Expected N
    N[g] ~ dpois(lambda.N[g]) #realized N in state space.
    
    #plug in shared df parameter for each session. Must use lam0[g] and sigma[g] here for custom update.
    #alternatively, can be estimated separately or with random effects.
    lam0[g] <- lam0.fixed
    sigma[g] <- sigma.fixed
    theta.thin[g] <- theta.thin.fixed
    
    #1) marked individuals
    for(i in 1:n.marked[g]){
      #uniform assumption here
      s.m[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s.m[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam.m[g,i,1:J[g]] <- GetDetectionRate(s=s.m[g,i,1:2],X=X[g,1:J[g],1:2],J=J[g],sigma=sigma[g],lam0=lam0[g],z=1)
      y.mID[g,i,1:J[g]] ~ dPoissonVector(lam.m[g,i,1:J[g]]*K1D[g,1:J[g]]*theta.thin[g],z=z.m[g,i]) #marked and identified detections
    }
    #marked with no ID detections
    bigLam.marked[g,1:J[g]] <- GetbigLam(lam=lam.m[g,1:n.marked[g],1:J[g]],z=z.m[g,1:n.marked[g]])
    lam.mnoID[g,1:J[g]] <- bigLam.marked[g,1:J[g]]*K1D[g,1:J[g]]*(1-theta.thin[g])
    y.mnoID[g,1:J[g]] ~ dPoissonVector(lam.mnoID[g,1:J[g]],z=1) #plug in z=1 to reuse dPoissonVector
    
    #2) all individuals treated as unmarked
    for(i in 1:M[g]){
      #dunif() here implies uniform distribution within a grid cell
      #also tells nimble s's are in continuous space, not discrete
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      #get cell s_i lives in using look-up table
      s.cell[g,i] <- cells[g,trunc(s[g,i,1]/res[g])+1,trunc(s[g,i,2]/res[g])+1]
      #categorical likelihood for this cell, equivalent to zero's trick
      #also disallowing s's in non-habitat
      dummy.data[g,i] ~ dCell(pi.cell[g,s.cell[g,i]])
      lam[g,i,1:J[g]] <- GetDetectionRate(s=s[g,i,1:2],X=X[g,1:J[g],1:2],J=J[g],sigma=sigma[g],lam0=lam0[g],z=z[g,i])
    }#custom Metropolis-Hastings update for N.M/z[1:n.marked] 
    bigLam.all[g,1:J[g]] <- GetbigLam(lam=lam[g,1:M[g],1:J[g]],z=z[g,1:M[g]])
    y.all[g,1:J[g]] ~ dPoissonVector(bigLam.all[g,1:J[g]]*K1D[g,1:J[g]],z=1) #plug in z=1 to reuse dPoissonVector
    
    #If you have telemetry
    for(i in 1:n.tel.inds[g]){
      for(m in 1:n.locs.ind[g,i]){
        locs[g,tel.inds[g,i],m,1] ~ dnorm(s.m[g,tel.inds[g,i],1],sd=sigma[g])
        locs[g,tel.inds[g,i],m,2] ~ dnorm(s.m[g,tel.inds[g,i],2],sd=sigma[g])
      }
    }
  }
})# end model