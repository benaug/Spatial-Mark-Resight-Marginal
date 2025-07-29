NimModel <- nimbleCode({
  #detection function priors - shared across sessions
  p0.fixed ~ dunif(0,1)
  lam0.fixed ~ dunif(0,15)
  sigma.fixed ~ dunif(0,10)
  
  #sample type observation probabilities, shared over sessions.
  alpha.marked[1] <- 1
  alpha.marked[2] <- 1
  alpha.marked[3] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  theta.marked.fixed[1:3] ~ ddirch(alpha.marked[1:3])
  theta.unmarked.fixed[1] <- 0
  theta.unmarked.fixed[2:3] ~ ddirch(alpha.unmarked[1:2])
  
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
    p0[g] <- p0.fixed
    lam0[g] <- lam0.fixed
    sigma[g] <- sigma.fixed
    
    #sample type observation model priors (Dirichlet)
    #If not shared across sessions, use this
    #into g indices here
    # alpha.marked[g,1] <- 1
    # alpha.marked[g,2] <- 1
    # alpha.marked[g,3] <- 1
    # alpha.unmarked[g,1] <- 1
    # alpha.unmarked[g,2] <- 1
    # theta.marked[g,1:3] ~ ddirch(alpha.marked[g,1:3])
    # theta.unmarked[g,1] <- 0
    # theta.unmarked[g,2:3] ~ ddirch(alpha.unmarked[g,1:2])
    #if shared over sessions use this
    theta.marked[g,1:3] <- theta.marked.fixed[1:3]
    theta.unmarked[g,1:3] <- theta.unmarked.fixed[1:3]
    
    #1 Marking Process
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
      pd[g,i,1:J.mark[g]] <- GetDetectionProb(s=s[g,i,1:2],X=X.mark[g,1:J.mark[g],1:2],J=J.mark[g],
                                              sigma=sigma[g],p0=p0[g],z=z[g,i])
      y.mark[g,i,1:J.mark[g]] ~ dBinomialVector(pd[g,i,1:J.mark[g]],K1D=K1D.mark[g,1:J.mark[g]],z=z[g,i])
    }
    
    #2 Sighting Process
    for(i in 1:M[g]){
      lam[g,i,1:J.sight[g]] <- GetDetectionRate(s=s[g,i,1:2],X=X.sight[g,1:J.sight[g],1:2],J=J.sight[g],
                                                sigma=sigma[g],lam0=lam0[g],z=z[g,i])
    }
    #2a Marked with ID detections
    for(i in 1:n.marked[g]){
      for(k in 1:K.sight[g]){
        #marked and identified detections. only count when marked.status[i,k] = 1
        y.mID[g,i,1:J.sight[g],k] ~ dPoissonVector(lam[g,i,1:J.sight[g]]*K2D.sight[g,1:J.sight[g],k]*
                                                     marked.status[g,i,k]*theta.marked[g,1],z=z[g,i])
      }
    }
    #Unidentified detections by type
    #2b marked with no ID detections
    for(k in 1:K.sight[g]){
      bigLam.marked[g,1:J.sight[g],k] <- GetbigLam(lam=lam[g,1:n.marked[g],1:J.sight[g]],z=z[g,1:n.marked[g]]*
                                                     marked.status[g,1:n.marked[g],k])
      lam.mnoID[g,1:J.sight[g],k] <- bigLam.marked[g,1:J.sight[g],k]*K2D.sight[g,1:J.sight[g],k]*theta.marked[g,2]
      y.mnoID[g,1:J.sight[g],k] ~ dPoissonVector(lam.mnoID[g,1:J.sight[g],k],z=1) #plug in z=1 to reuse dPoissonVector
    }
    
    #2c unmarked detections
    bigLam.unmarked[g,1:J.sight[g]] <- GetbigLam(lam=lam[g,(n.marked[g]+1):M[g],1:J.sight[g]],z=z[g,(n.marked[g]+1):M[g]])
    for(k in 1:K.sight[g]){
      bigLam.premarked[g,1:J.sight[g],k] <- GetbigLam(lam=lam[g,1:n.marked[g],1:J.sight[g]],
                                                      z=z[g,1:n.marked[g]]*(1-marked.status[g,1:n.marked[g],k]))
      lam.um[g,1:J.sight[g],k] <- bigLam.unmarked[g,1:J.sight[g]]*K2D.sight[g,1:J.sight[g],k]*theta.unmarked[g,2] + #always unmarked individuals
        bigLam.premarked[g,1:J.sight[g],k]*K2D.sight[g,1:J.sight[g],k]*theta.unmarked[g,2] #marked inds when not marked
      y.um[g,1:J.sight[g],k] ~ dPoissonVector(lam.um[g,1:J.sight[g],k],z=1) #plug in z=1 to reuse dPoissonVector
    }
    
    #2d unknown marked status
    for(k in 1:K.sight[g]){
      lam.unk[g,1:J.sight[g],k] <- bigLam.marked[g,1:J.sight[g],k]*K2D.sight[g,1:J.sight[g],k]*theta.marked[g,3] + #marked individuals when marked
        bigLam.premarked[g,1:J.sight[g],k]*K2D.sight[g,1:J.sight[g],k]*theta.unmarked[g,3] + #marked individuals when not marked
        bigLam.unmarked[g,1:J.sight[g]]*K2D.sight[g,1:J.sight[g],k]*theta.unmarked[g,3] #always unmarked individuals
      y.unk[g,1:J.sight[g],k] ~ dPoissonVector(lam.unk[g,1:J.sight[g],k],z=1) #plug in z=1 to reuse dPoissonVector
    }
    
    #If you have telemetry
    for(i in 1:n.tel.inds[g]){
      for(m in 1:n.locs.ind[g,i]){
        locs[g,i,m,1] ~ dnorm(s[g,tel.inds[g,i],1],sd=sigma[g])
        locs[g,i,m,2] ~ dnorm(s[g,tel.inds[g,i],2],sd=sigma[g])
      }
    }
  }
})# end model