library(coda)
burnin <- 500
n.datasets <- 100
N.session <- 3
true <- c(rep(0.5,N.session),rep(exp(-0.5),N.session),
          rep(NA,N.session),0.5,rep(NA,N.session),
          0.5,0.75,0.15,0.1,0,0.75,0.25)
n.par <- length(true)
true <- matrix(true,nrow=n.datasets,ncol=n.par,byrow=TRUE)
ests1 <- ests2 <- ests3 <- bias1 <- bias2 <- bias3 <- cover <- sds <- matrix(NA,n.datasets,n.par)
N.max <- time <- rep(NA,n.datasets)
HPDs <- array(NA,c(n.datasets,n.par,2))
for(d in 1:n.datasets){
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Multisession_tel_S1_datasets")
  filename <- paste0("S1_data_",d,".RData")
  load(filename)
  setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing/Multisession_tel_S1_fits")
  filename <- paste0("S1_",d,".chain_1.RData")
  load(filename)
  post <- mcmc(out$mvSamples[-c(1:burnin),])
  # plot(post)
  time[d] <- out$time
  true[d,7:9] <- sapply(data,function(x){x$N})
  true[d,11:13] <- sapply(data,function(x){x$lambda.N})
  N.max[d] <- max(post[,6:9])
  ests1[d,] <- colMeans(post)
  ests2[d,] <- MCMCglmm::posterior.mode(post)
  ests3[d,] <- apply(post,2,quantile,probs=0.5)
  #N and lambda need higher bandwidths
  # ests2[d,c(7:9)] <- MCMCglmm::posterior.mode(post[,c(7:9)],adjust=1)
  bias1[d,] <- (ests1[d,] - true[d,])/true[d,]
  bias2[d,] <- (ests2[d,] - true[d,])/true[d,]
  bias3[d,] <- (ests3[d,] - true[d,])/true[d,]
  HPDs[d,,] <- HPDinterval(post)
  sds[d,] <- apply(post,2,sd)
  cover[d,] <- 1*(HPDs[d,,1]<=true[d,]&HPDs[d,,2]>=true[d,])
}


par(mfrow=c(1,1),ask=FALSE)
par.names <- colnames(post)

colMeans(true)
colMeans(ests1) #means
colMeans(ests2,na.rm=TRUE) #modes
colMeans(cover,na.rm=TRUE)
100*colMeans(bias1)
100*colMeans(bias2)
#CVs
CVs <- 100*colMeans(sds/abs(ests2))

#modes
summary2 <- rbind(colMeans(true),colMeans(ests2),100*colMeans(bias2),colMeans(cover),CVs)
colnames(summary2) <- par.names
rownames(summary2) <- c("true","est","rel.bias","cover","CV")
round(summary2,3)


sort(N.max)
mean(time)
