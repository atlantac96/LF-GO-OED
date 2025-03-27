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

#Generate 1000 samples from prior of theta
set.seed(123)
epsilon<-rnorm(1000,0,0.01)
theta<-runif(1000,0,1)
#Goal= parameter inference
z<-theta

#discrete design variable
design<-seq(0,1,length.out=21)

#samples of observed data from prior samples
y<-matrix(NA, nrow=length(design), ncol=1000)
for (k in 1: length(design)){
  d<- design[k]
  y[k,]<- theta^3*d^2+theta* exp(-abs(0.2-d))+epsilon
}

#To store the posterior samples of theta
theta.abc<-list()
zabc<-list()

rat<-matrix(NA, nrow=length(design), ncol=length(z))
for(k in 1:21){
    cat(" Starting Replicate with design point",k,
        "at", print(Sys.time()), "\n")
  #store posterior samples for each discrete design value
    theta.abc[[k]]<- list() 
    #Run the posterior analysis and density ratio in parallel
    library(foreach)
    library(doParallel)
    rat[k,]<-foreach(i = 1:length(y[k,]), .combine=cbind) %dopar% {
      #do ABC posterior
      library(abc)
      # do cross validation to find tolerance rate,  expensive step
      # roughly 6 mins for 1000 nval
      cvtol<-cv4abc(param=theta, sumstat=y[k,], nval=1000,tols=c(0.1,0.3,0.5), method="neuralnet")
      #check difference between estim and true, select tolerance rate based on which is lower
      
      # OR use fixed tolerance after doing cross validation, after running above for one design value
      #3 seconds for ABC
      check<- abc(y[k,i],theta, y[k,], tol=0.3, method="neuralnet", sizenet=120)

      theta.abc[[k]][[i]]<-check$adj.values # or set it to cvtol$estim
      zabc[[k]]<- check$adj.values # since goal=paramter inference
      
      #generate prior samples for denominator of utility evaluation
      set.seed(i)
      thetai<-runif(nrow(check$adj.values),0,1)
      zi<-thetai 
      
      #computing the density ratio
      library(densratio)
      densratio_obj <- densratio::densratio(check$adj.values, zi,method="RuLSIF", 
                                            alpha=0.0,kernel_num=300, fold=10)
      #use options below to run Rulsif and hyperparameter validation 
      #alpha=0.2 ,sigma = 10 ^ seq(-5, 1, length.out = 50),
      #lambda = 10 ^ seq(-5, 1, length.out = 50))
      valuea<-densratio_obj$compute_density_ratio(z)[i]
      
    }
}

ratna<- log(rat)
ratna[!is.finite(ratna)] <- NA
ut<-rowMeans(ratna, na.rm=TRUE)
plot(design,ut, type="b")

#analytical evaluation
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
points(design, ut.ana, col="red", type="b")

#LFIRE
set.seed(123)
epsilon<-rnorm(1000,0,0.01)
theta<-runif(1000,0,1)
z<-theta

design<-seq(0,1,length.out=21)
y<-matrix(NA, nrow=length(design), ncol=1000)
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
    eplik<- rnorm(1000,0,0.01)
    ylik<-theta[i]^3*design[k]^2+theta[i]* exp(-abs(0.2-design[k]))+eplik
    
    set.seed(2*i)
    thetam<-runif(1000,0,1)
    epmarg<- rnorm(1000,0,0.01)
    ymarg<-thetam^3*design[k]^2+thetam* exp(-abs(0.2-design[k]))+epmarg
    
    t_t <- rep(1, 1000)
    t_m <- rep(0, 1000)
    
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


