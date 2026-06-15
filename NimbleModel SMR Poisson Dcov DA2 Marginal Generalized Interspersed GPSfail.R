NimModel <- nimbleCode({
  #baseline density
  D0 ~ dunif(0,100) #uninformative, diffuse dnorm on log scale can cause neg bias
  # D.beta0 ~ dnorm(0,sd=10)
  # D0 <- exp(D.beta0)
  D.beta1 ~ dnorm(0,sd=10)

  #detection function priors
  p0 ~ dunif(0,1) #marking process
  lam0 ~ dunif(0,15) #sighting proces
  sigma ~ dunif(0,10) #marking and sighting process
  #sample type observation model priors (Dirichlet)
  alpha.marked1[1] <- 1
  alpha.marked1[2] <- 1
  alpha.marked1[3] <- 1
  alpha.marked2[1] <- 1
  alpha.marked2[2] <- 1
  alpha.unmarked[1] <- 1
  alpha.unmarked[2] <- 1
  # mark type 1: working GPS collar: all obs types possible
  theta.marked1[1:3] ~ ddirch(alpha.marked1[1:3])
  # mark type 2: dead GPS collar: cannot produce individual ID
  theta.marked2[1] <- 0
  theta.marked2[2:3] ~ ddirch(alpha.marked2[1:2])
  # truly unmarked / not-yet-marked: cannot produce individual ID
  theta.unmarked[1] <- 0
  theta.unmarked[2:3] ~ ddirch(alpha.unmarked[1:2])

  #Density model
  D.intercept <- D0*cellArea
  lambda.cell[1:n.cells] <- InSS[1:n.cells]*exp(D.beta1*D.cov[1:n.cells])
  pi.cell[1:n.cells] <- lambda.cell[1:n.cells]/pi.denom #expected proportion of total N in cell c
  pi.denom <- sum(lambda.cell[1:n.cells])
  lambda.N <- D.intercept*pi.denom #Expected N
  N ~ dpois(lambda.N) #realized N in state space.
  
  #1 Marking Process
  for(i in 1:M){
    s[i,1] ~ dunif(xlim[1],xlim[2])
    s[i,2] ~ dunif(ylim[1],ylim[2])
    #get cell s_i lives in using look-up table
    s.cell[i] <- cells[trunc(s[i,1]/res)+1,trunc(s[i,2]/res)+1]
    #categorical likelihood for this cell, equivalent to zero's trick
    #also disallowing s's in non-habitat
    dummy.data[i] ~ dCell(pi.cell[s.cell[i]])
    pd[i,1:J.mark] <- GetDetectionProb(s=s[i,1:2],X=X.mark[1:J.mark,1:2],J=J.mark,sigma=sigma,p0=p0,z=z[i])
    y.mark[i,1:J.mark] ~ dBinomialVector(pd[i,1:J.mark],K1D=K1D.mark[1:J.mark],z=z[i])
  }
  
  #2 Sighting Process
  for(i in 1:M){
    lam[i,1:J.sight] <- GetDetectionRate(s=s[i,1:2],X=X.sight[1:J.sight,1:2],J=J.sight,sigma=sigma,lam0=lam0,z=z[i])
  }
  #2a Marked with ID detections
  for(i in 1:n.marked){
    for(k in 1:K.sight){
      lam.mID[i,1:J.sight,k] <- lam[i,1:J.sight]*K2D.sight[1:J.sight,k]*
        (mark.class1[i,k]*theta.marked1[1] + mark.class2[i,k]*theta.marked2[1])
      y.mID[i,1:J.sight,k] ~ dPoissonVector(lam.mID[i,1:J.sight,k],z=z[i])
    }
  }
  #2b marked with no ID detections
  for(k in 1:K.sight){
    bigLam.marked1[1:J.sight,k] <- GetbigLam(lam=lam[1:M,1:J.sight],
                                             z=z[1:M]*mark.class1[1:M,k])
    bigLam.marked2[1:J.sight,k] <- GetbigLam(lam=lam[1:M,1:J.sight],
                                             z=z[1:M]*mark.class2[1:M,k])
    lam.mnoID[1:J.sight,k] <- K2D.sight[1:J.sight,k]*
      (bigLam.marked1[1:J.sight,k]*theta.marked1[2] +
         bigLam.marked2[1:J.sight,k]*theta.marked2[2])
    y.mnoID[1:J.sight,k] ~ dPoissonVector(lam.mnoID[1:J.sight,k],z=1)
  }
  #2c unmarked detections
  for(k in 1:K.sight){
    bigLam.unmarked[1:J.sight,k] <- GetbigLam(lam=lam[1:M,1:J.sight],
                                              z=z[1:M]*mark.class0[1:M,k])
    lam.um[1:J.sight,k] <- bigLam.unmarked[1:J.sight,k]*
      K2D.sight[1:J.sight,k]*theta.unmarked[2]
    y.um[1:J.sight,k] ~ dPoissonVector(lam.um[1:J.sight,k],z=1)
  }
  #2d unknown marked status
  for(k in 1:K.sight){
    lam.unk[1:J.sight,k] <- K2D.sight[1:J.sight,k]*
      (bigLam.marked1[1:J.sight,k]*theta.marked1[3] +
         bigLam.marked2[1:J.sight,k]*theta.marked2[3] +
         bigLam.unmarked[1:J.sight,k]*theta.unmarked[3])
    y.unk[1:J.sight,k] ~ dPoissonVector(lam.unk[1:J.sight,k],z=1)
  }
  
  #If you have telemetry
  for(i in 1:n.tel.inds){
    for(m in 1:n.locs.ind[i]){
      locs[i,m,1] ~ dnorm(s[tel.inds[i],1],sd=sigma)
      locs[i,m,2] ~ dnorm(s[tel.inds[i],2],sd=sigma)
    }
  }
})# end model