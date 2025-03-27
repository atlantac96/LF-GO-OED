rm(list=ls())
library(curl)
library(dplyr)
library(abc)
library(densratio)
library(doParallel)
totalCores = detectCores()
cluster <- makeCluster(totalCores[1]) 
registerDoParallel(cluster)
library(ParBayesianOptimization)
library(foreach)


start<-Sys.time()
ndim<-15
param_space <- list()
for (i in 1:ndim) {
  param_space[[paste0("params", i)]] <- c(0,1)
}
objective_function <- function(... ) {
  param_names <- names(list(...))
  param_values <- unlist(list(...))
  ndim<-20
  set.seed(123)
  epsilon<-matrix(rnorm(500*ndim,0,0.01), nrow=500, ncol=ndim)
  theta<-matrix(runif(500*ndim,0,1), nrow=500,ncol=ndim)
  calculate_z <- function(theta_row) {
    sum((1 - theta_row[-length(theta_row)])^2 + 5 * (theta_row[-1] - theta_row[-length(theta_row)])^2)
  }
  z <- apply(theta, 1, calculate_z)
  
  y<-list()
  theta.abc<-list()
  zabc<-list()
  
  rat<-rati<-rep(NA,length(z)) #matrix(NA, nrow=nrow(design), ncol=length(z))
  
  
    y<-matrix(NA, nrow=500, ncol=ndim)
    for (i in 1:ndim) {
      y[,i] <- theta[,i]^3 * param_values[i]^2 + 
        sum(theta[, -i] * exp(-abs(0.2 - param_values[-i]))) + epsilon[,i]
    }
    
  library(foreach)
  library(doParallel)
  rat<-foreach(i = 1:nrow(y), .combine=cbind) %dopar% {
    calculate_z <- function(theta_row) {
      sum((1 - theta_row[-length(theta_row)])^2 + 5 * (theta_row[-1] - theta_row[-length(theta_row)])^2)
    }
    library(abc)
    check<- abc(y[i,],theta, y, tol=0.1, method="neuralnet")
    theta.abc<-check$adj.values
    zabc<- apply(check$adj.values, 1, calculate_z)
    
    set.seed(i)
    thetai<-matrix(runif(ndim*nrow(check$adj.values),0,1), 
                   nrow=nrow(check$adj.values),ncol=ndim)
    zi<-apply(check$adj.values, 1, calculate_z)
    library(densratio)
    densratio_obj <- densratio::densratio(zabc, zi)#sigma=0.01,lambda=0.01
    rati[i]<-densratio_obj$compute_density_ratio(z)[i]
  }
  ratna<- log(rat)
  ratna[!is.finite(ratna)] <- NA
  ut <- rowMeans(ratna, na.rm = TRUE)
  return(list(Score = ut))
  
}
clusterExport(cluster, varlist=c("objective_function"), envir=.GlobalEnv)
set.seed(42)
res<- bayesOpt(FUN=objective_function, bounds=param_space,
                    plotProgress = TRUE,
                    initPoints =50,
                    iters.n=25,
                    acq = "ei", 
                    nugget=0.01,
                    #convThresh = 90,
                    parallel = TRUE, verbose=2 , 
                    otherHalting = list(timeLimit = Inf, minUtility=-0.001) )
end<- Sys.time()
timetaken<- end-start
getLocalOptimums(res)
df<-res$scoreSummary
maxut[k]<-max(df$Score)
min(df$gpUtility, na.rm=TRUE)
getBestPars(res)
maxut
# 2 dim, takes 32.22053 mins with fixed alpha and lanbda at 0.01, 
# 2 dim takes   mins with densratio cross validation 1.092907 hours
#3 dim  of 1.077225 hours
#4 dim Time difference of 1.115706 hours

library(ggplot2)
pdf("dimvstime.pdf", height=4, width=8)
par(mar=c(5,6,4,1)+.1)
plot(data$Dimension, data$TimeTaken, type = "b",cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1,
     xlab = "Dimension",
     ylab = "Time taken (hrs)") 
dev.off()
