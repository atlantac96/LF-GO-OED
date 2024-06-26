rm(list=ls())
library(curl)
library(dplyr)
library(abc)
library(densratio)
library(doParallel)
totalCores = detectCores()
cluster <- makeCluster(totalCores[1]) 
registerDoParallel(cluster)
library(foreach)

set.seed(123)
epsilon<-matrix(rnorm(20000,0,0.01), nrow=10000, ncol=2)
theta<-matrix(runif(20000,0,1), nrow=10000,ncol=2)
z<-(1-theta[,1])^2+(5*(theta[,2]-theta[,1]^2)^2) #rosenbrock

y<-list()
theta.abc<-list()
zabc<-list()

design<-expand.grid(seq(0,1,length.out=15),seq(0,1,length.out=15))
rat<-matrix(NA, nrow=nrow(design), ncol=length(z))
rati<-matrix(NA, nrow=nrow(design), ncol=length(z))
for(k in 1:nrow(design)){
  cat(" Starting Replicatewith design point",k,
      "at", print(Sys.time()), "\n")
  y[[k]]<- cbind(theta[,1]^3*design[k,1]^2+theta[,2]* exp(-abs(0.2-design[k,2]))+epsilon[,1],
                 theta[,2]^3*design[k,1]^2+theta[,1]* exp(-abs(0.2-design[k,2]))+epsilon[,2])
  
  theta.abc[[k]]<- list() 
  library(foreach)
  library(doParallel)
  rat[k,]<-foreach(i = 1:nrow(y[[k]]), .combine=cbind) %dopar% {
    library(abc)
    check<- abc(y[[k]][i,],theta, y[[k]], tol=0.01, method="neuralnet")
    theta.abc[[k]][[i]]<-check$adj.values
    zabc[[k]]<- (1-check$adj.values[,1])^2+(5*(check$adj.values[,2]-
                                            check$adj.values[,1]^2)^2)
    
    set.seed(i)
    thetai<-matrix(runif(2*nrow(check$adj.values),0,1), nrow=nrow(check$adj.values),ncol=2)
    zi<-(1-thetai[,1])^2+(5*(thetai[,2]-thetai[,1]^2)^2)
    library(densratio)
    densratio_obj <- densratio::densratio(zabc[[k]], zi)
                                       # alpha=0.0, kernel_num = 50)
                #sigma = 10 ^ seq(-5, 1, length.out = 100),
                 #lambda = 10 ^ seq(-5, 1, length.out = 100))
    rati[k,i]<-densratio_obj$compute_density_ratio(z)[i]
  }
}
ratna<- log(rat)
ratna[!is.finite(ratna)] <- NA
save(ut, file="2d2thetarosenbrock.RDa")
load("2d2thetarosenbrock.RDa") #this is the utility
design<-expand.grid(seq(0,1,length.out=15),seq(0,1,length.out=15))
dunique<-unique(design[,1])
utmatrix<-matrix(NA, nrow=length(dunique), ncol=length(dunique))
for( i in 1:length(ut)){
  utmatrix[which(dunique==design[i,1]),which(dunique==design[i,2])]<-ut[i]-0.6
}
cols = rev(colorRampPalette(c('blue', 'green'))(24))
pdf(file = "utility2d2theta.pdf",
    width = 8, 
    height = 8)
zlim=seq(1.3,3.8, by=0.8)
par(mar=c(5,6,4,1)+.1)

filled.contour(x=dunique[-1], y=dunique[-1], z=utmatrix[-1,-1], 
               zlim=zlim,
               plot.title = title(
                                  xlab=bquote(paste(italic(d[1]))), 
                                  ylab=bquote(paste(italic(d[2]))),
                                  cex.lab=2,
                                  cex.main=2),
               col=cols,
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(x=dunique[-1], y=dunique[-1], z=utmatrix[-1,-1], 
                         add = TRUE, labcex=1.5, nlevels=10
                 )
               }
)
dev.off()
#Stop cluster
stopCluster(cluster)
dnmcrosenbrock<-read.delim2("KL2D_funcPred_RR_ADBW", sep=",")
dnmcrosenbrock[ , 2:4] <- apply(dnmcrosenbrock[ , 2:4], 2, function(x) as.numeric(as.character(x)))
design<-expand.grid(seq(0,1,length.out=15),seq(0,1,length.out=15))
dunique<-unique(design[,1])
utmatrixdnmc<-matrix(NA, nrow=length(dunique), ncol=length(dunique))
for( i in 1:nrow(dnmcrosenbrock)){
  utmatrixdnmc[which(dunique==dnmcrosenbrock[i,3]),which(dunique==dnmcrosenbrock[i,2])]<-dnmcrosenbrock[i,4]
}
cols = rev(colorRampPalette(c('blue', 'green'))(24))
pdf(file = "utility2d2thetarosenbrockdnmc.pdf",
    width = 8, 
    height = 8)

par(mar=c(5,6,4,1)+.1)
filled.contour(x=dunique[-1], y=dunique[-1], z=utmatrixdnmc[-1,-1], 
               zlim=zlim,
               plot.title = title(
                                xlab=bquote(paste(italic(d[1]))), 
                                  ylab=bquote(paste(italic(d[2]))),
                                 cex.lab=2,
                                  cex.main=2),
               col=cols,
               plot.axes = {
                 axis(1)
                 axis(2)
                 contour(x=dunique[-1], y=dunique[-1], z=utmatrixdnmc[-1,-1], 
                         add = TRUE, nlevels=10, labcex = 1.5
                         )
               }
)
dev.off()



