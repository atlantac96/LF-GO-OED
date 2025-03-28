rm(list=ls())
library(curl)
library(dplyr)
library(abc)
library(densratio)
library(glmnet)
library(doParallel)
totalCores <- detectCores()
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)

start<-Sys.time()
set.seed(123)
epsilon<-rnorm(5000,0,0.01)
theta<-runif(5000,0,1)
z<-theta

design<-seq(0,1,length.out=21)
y<-matrix(NA, nrow=length(design), ncol=5000)
for (k in 1: length(design)){
  d<- design[k]
  y[k,]<- theta^3*d^2+theta* exp(-abs(0.2-d))+epsilon
}

theta.abc<-list()
zabc<-list()

rat<-matrix(NA, nrow=length(design), ncol=length(z))
for(k in 1:21){
  cat(" Starting Replicatewith design point",k,
      "at", print(Sys.time()), "\n")
  theta.abc[[k]]<- list() 
  library(foreach)
  library(doParallel)
  start<-Sys.time()
  rat[k,]<-foreach(i = 1:length(y[k,]), .combine=cbind) %dopar% {
    library(abc)
    check<- abc(y[k,i],theta, y[k,], tol=0.3, method="neuralnet")
    theta.abc[[k]][[i]]<-check$adj.values
    zabc[[k]]<- check$adj.values
    set.seed(i)
    thetai<-runif(nrow(check$adj.values),0,1)
    zi<-thetai 
    library(densratio)
    
    densratio_obj <- densratio::densratio(check$adj.values, zi)
    #alpha=0.2 ,sigma = 10 ^ seq(-5, 1, length.out = 50),
    #lambda = 10 ^ seq(-5, 1, length.out = 50))
    valuea<-densratio_obj$compute_density_ratio(z)[i]
    
  }
  stop<-Sys.time()
  timetaken<-stop-start
}
ratna<- log(rat)
ratna[!is.finite(ratna)] <- NA
ut<-rowMeans(ratna, na.rm=TRUE)
end<-Sys.time()
timetaken<-end-start
plot(design,ut, type="b")


#analytcal evaluation
start<-Sys.time()
set.seed(123)
epsilon<-rnorm(1000,0,0.01)
theta<-runif(1000,0,1)
z<-theta #sin(theta)+theta*exp(theta+abs(0.5-theta))

design<-seq(0,1,length.out=21)
y<-matrix(NA, nrow=length(design), ncol=1000)
xtest<-c()
for (k in 1: length(design)){
  d<- design[k]
  y[k,]<- theta^3*d^2+theta* exp(-abs(0.2-d))+epsilon
  xtest<-rbind(xtest,cbind(rep(d,1000),y[k,],z))
}
eval<-function(d,yv,thetav){
  py.thetad<-(1/(sqrt(2*pi*10^{-4})))* exp(-(yv-thetav^3*d^2-thetav*exp(-abs(0.2-d)))^2/(2*10^{-4}))
  return(py.thetad)
}
ut.ana<-c()
colnames(xtest)[1:3]=c("d", "yv", "thetav")
lik<-temp<-c()
for(k in 1:length(design)){
  ff<-xtest[xtest[,1]==design[k],]
  for(j in 1:nrow(ff)){
    lik[j]<-eval(ff[j,1],ff[j,2],ff[j,3])
    marg<-eval(ff[j,1],ff[j,2],theta)
    temp[j]<-lik[j]/mean(marg)
  }
  ut.ana[k]<-mean(log(temp))
}
plot(design, ut.ana, col="red", type="b")
end<-Sys.time()
timetaken<-end-start

#lfire 
set.seed(123)
epsilon<-rnorm(5000,0,0.01)
theta<-runif(5000,0,1)
z<-theta

design<-seq(0,1,length.out=21)
y<-matrix(NA, nrow=length(design), ncol=5000)
for (k in 1: length(design)){
  d<- design[k]
  y[k,]<- theta^3*d^2+theta* exp(-abs(0.2-d))+epsilon
}


rat<-matrix(NA, nrow=length(design), ncol=length(z))
for(k in 1:21){
  cat(" Starting Replicate with design point",k,
      "at", print(Sys.time()), "\n")
  theta.abc[[k]]<- list() 
  library(foreach)
  library(doParallel)
  rat[k,]<-foreach(i = 1:length(y[k,]), .combine=cbind) %dopar% {
    set.seed(i)
    eplik<- rnorm(5000,0,0.01)
    ylik<-theta[i]^3*design[k]^2+theta[i]* exp(-abs(0.2-design[k]))+eplik
    
    set.seed(2*i)
    thetam<-runif(5000,0,1)
    epmarg<- rnorm(5000,0,0.01)
    ymarg<-thetam^3*design[k]^2+thetam* exp(-abs(0.2-design[k]))+epmarg
    
    t_t <- rep(1, 5000)
    t_m <- rep(0, 5000)
    
    Y <-as.data.frame(c(ylik, ymarg))
    colnames(Y) <- "y"
    class <- as.factor(c(t_t, t_m))
    
    model <- glm(class ~ ., data = Y, family = binomial())
    print(summary(model)$coefficients)
    
    predicted <- ifelse(predict(model, type = "response") > 0.5, 1, 0)
    confusion_matrix <- table(Actual =class, Predicted = predicted)
    print(confusion_matrix)
    
    ydf<-as.data.frame(y[k,])
    colnames(ydf) <- "y"
    probabilities <- predict(model, newdata =ydf, 
                             type = "response")[i]
  }
}
ratna<-log(rat/(1-rat))
ratna[!is.finite(ratna)] <- NA
utlfire<-rowMeans(ratna, na.rm=TRUE)
plot(design,utlfire)


library(ggplot2)
pdf("toy1d.pdf", height=4, width=8)
par(mar=c(5,6,4,1)+.1)
plot(design, ut, type = "b",cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1,
     xlab = expression(paste(italic("d"))),
     ylab = bquote(paste(italic(U[z](d))))) 
points(design, ut.ana, col="red", type="b")
legend("topleft", legend = c("LF-GO-OED", "Nested MC"),
       col = c("black", "red"), lty =1, cex = 1)
dev.off()


