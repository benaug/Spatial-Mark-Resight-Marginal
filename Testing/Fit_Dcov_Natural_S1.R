setwd("D:/Sync/Cornell/Github/Spatial-Mark-Resight-Marginal/Testing")
dir.create("Dcov_Natural_S1_fits")
n.datasets <- 100
n.cores <- 18

seeds <- round(runif(n.datasets,1,10000))

library(snow)
library(doSNOW)
library(foreach)
cl.tmp  <-  makeCluster(rep("localhost",n.cores), type="SOCK")
registerDoSNOW(cl.tmp)
out <- foreach(d=1:n.datasets) %dopar% {
  system2("C:\\Users\\ba378\\AppData\\Local\\Programs\\R\\R-4.4.0\\bin\\Rscript.exe", 
          c("D:\\Sync\\Cornell\\Github\\Spatial-Mark-Resight-Marginal\\Testing\\Fit_Dcov_Natural_S1_single.R",d=d,seeds[d]), invisible=TRUE)
}
stopCluster(cl.tmp)
