rm(list=ls())
library(deSolve)
library(doParallel)
totalCores <- detectCores()
cluster <- makeCluster(totalCores[1]) 
registerDoParallel(cluster)
library(foreach)
library(truncnorm)
library(odin)

library(jsonlite)

data <- fromJSON("data.json")
design<- data$ts
num_replicates <- 1
ut_matrix <- matrix(NA, nrow = length(design), ncol = num_replicates)



for(replicate_num in 1:num_replicates){
  
  
  library(doParallel)
  cat("Replicate:", replicate_num, "\n")
  
  
  # Priors
  beta = data$prior_samples[,1]
  gamma = data$prior_samples[,2]
  parameters_values <- priorabc<-cbind(beta, gamma)
  
  z_recov_prior <- parameters_values
  sir_results<-data$ys
  
  
  trans <- vector("list", length = length(design))
  result_matrix <- matrix(NA, nrow = nrow(priorabc), ncol = 2)
  
  for (s in 1: length(design)) {
    for (i in 1:nrow(priorabc)){
      result_matrix[i, 1] <-sir_results[s,i,1]
      result_matrix[i, 2] <-sir_results[s,i,2]
    }
    trans[[s]] <- result_matrix
  }
  
  
  z.abc <- list()
  rat <- rati <- matrix(NA, nrow = length(design), ncol = nrow(parameters_values))
  for (k in 2:length(design)) {
    cat(" Starting Replicate:", replicate_num, "with design point",k,
        "at", print(Sys.time()), "\n")
    z.abc[[k]] <- list()
    rat[k, ] <- foreach(i = 1:nrow(parameters_values), .combine = cbind) %dopar% {
      
      library(deSolve)
      library(truncnorm)
      library(abc)
      check<-abc(trans[[k]][i,],priorabc, trans[[k]],
                 method="neuralnet", tol=0.01)
      summary(check$unadj.values)
      priorabc[i,]
      
      z.abc[[k]][[i]] <- check$adj.values[complete.cases(check$adj.values), ]
      
      set.seed(i)
      thetai <- matrix(NA, nrow = nrow(check$adj.values), ncol = 2)
      thetai[, 1] =  rlnorm(nrow(check$adj.values), meanlog = 0.50, sdlog = 0.50)
      thetai[, 2] = rlnorm(nrow(check$adj.values), meanlog = 0.10, sdlog = 0.50)
      colnames(thetai) <- c("beta", "gamma")
      zi <- thetai[complete.cases(thetai), ]
      
      library(densratio)
      set.seed(56)
      densratio_obj <- densratio::densratio(z.abc[[k]][[i]], zi, 
                                          method = "RuLSIF",
                                           alpha = 1,
                                          sigma = 10 ^ seq(-5, 3, length.out = 10),
                                          lambda = 10 ^ seq(-5, 3, length.out = 10))
      rati[k, i] <- densratio_obj$compute_density_ratio(z_recov_prior)[i]
      # }
    }
  }
  ratna <- log(rat)
  ratna[!is.finite(ratna)] <- NA
  
  ut <- rowMeans(ratna, na.rm = TRUE)
  ut[k]
  ut_matrix[, replicate_num] <- ut
}
library(ggplot2)
pdf("appendixsir.pdf", height=4, width=8)
plot(design, ut, type = "b",cex.lab=1.2, cex.axis=1, cex.main=1, cex.sub=1,
     xlab = expression(paste("Design, ", italic("d"))),
     ylab = bquote(paste("Utility, ", italic(U[z](d))))) 
dev.off()

plot(design,ut)
save(ut_matrix, file="utmatrixstochasticsirparaminf.RDa")

stopCluster(cluster)
library(ggplot2)
load("utmatrixstochasticsirparaminf.RDa")
mean_row <- rowMeans(ut_matrix, na.rm=TRUE)
sd_row <- apply(ut_matrix, 1, sd, na.rm=TRUE)
df <- data.frame(sequence = 1:nrow(ut_matrix), mean = mean_row, sd = sd_row)
pdf("sirparaminf.pdf", height=4, width=8)
ggplot(data=df, aes(x=design, y=mean)) + 
  labs(x = expression(paste("Design, ", italic("d"))), 
       y = bquote(paste("Utility, ", italic(U[z](d)))), 
       color = "") +  theme_minimal() +
  geom_line() + theme(text=element_text(size=20))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha = 0.2)
dev.off()
design[which.max(df$mean)] #0.4
df$mean[which.max(df$mean)] #2.066
design[which.min(df$sd)]
df$sd[which.min(df$sd)]
df$sd[which.max(df$sd)]
ut
[1]      NaN 1.064661 1.083892 1.101734 1.109695 1.119411 1.119551 1.150567
[9] 1.151369 1.144831 1.171915 1.185702 1.194374 1.194145 1.199809 1.225106
[17] 1.228116 1.229630 1.239621 1.239126 1.265188 1.267206 1.287114 1.287287
[25] 1.294814 1.287985 1.308489 1.318137 1.313478 1.328694

#lfire:
data <- fromJSON("data.json")
design<- data$ts
num_replicates <- 1
ut_matrix <- matrix(NA, nrow = length(design), ncol = num_replicates)



for(replicate_num in 1:num_replicates){
  
  
  library(doParallel)
  cat("Replicate:", replicate_num, "\n")
  
  
  # Priors
  beta = data$prior_samples[,1]
  gamma = data$prior_samples[,2]
  parameters_values <- priorabc<-cbind(beta, gamma)
  
  z_recov_prior <- parameters_values
  sir_results<-data$ys
  
  
  trans <- vector("list", length = length(design))
  result_matrix <- matrix(NA, nrow = nrow(priorabc), ncol = 2)
  
  for (s in 1: length(design)) {
    for (i in 1:nrow(priorabc)){
      result_matrix[i, 1] <-sir_results[s,i,1]
      result_matrix[i, 2] <-sir_results[s,i,2]
    }
    trans[[s]] <- result_matrix
  }
  
  
  z.abc <- list()
  rat <- rati <- matrix(NA, nrow = length(design), ncol = nrow(parameters_values))
  for (k in 2:length(design)) {
    cat(" Starting Replicate:", replicate_num, "with design point",k,
        "at", print(Sys.time()), "\n")
    z.abc[[k]] <- list()
    rat[k, ] <- foreach(i = 1:nrow(parameters_values), .combine = cbind) %dopar% {
      
      library(deSolve)
      library(truncnorm)
      library(abc)
      check<-abc(trans[[k]][i,],priorabc, trans[[k]],
                 method="neuralnet", tol=0.01)
      summary(check$unadj.values)
      priorabc[i,]
      
      z.abc[[k]][[i]] <- check$adj.values[complete.cases(check$adj.values), ]
      
      set.seed(i)
      thetai <- matrix(NA, nrow = nrow(check$adj.values), ncol = 2)
      thetai[, 1] =  rlnorm(nrow(check$adj.values), meanlog = 0.50, sdlog = 0.50)
      thetai[, 2] = rlnorm(nrow(check$adj.values), meanlog = 0.10, sdlog = 0.50)
      colnames(thetai) <- c("beta", "gamma")
      zi <- thetai[complete.cases(thetai), ]
      t_t <- rep(1, nrow(check$adj.values))
      t_m <- rep(0, nrow(check$adj.values))
      
      Y <-as.data.frame(rbind(z.abc[[k]][[i]], thetai))
      #colnames(Y) <- c("y")
      class <- as.factor(c(t_t, t_m))
      
      model <- glm(class ~ ., data = Y, family = binomial())
      print(summary(model)$coefficients)
      
      predicted <- ifelse(predict(model, type = "response") > 0.5, 1, 0)
      confusion_matrix <- table(Actual =class, Predicted = predicted)
      print(confusion_matrix)
      
      ydf<-as.data.frame(priorabc)
      probabilities <- predict(model, newdata =ydf, 
                               type = "response")[i]
    }
  }
  ratna <- log(rat)
  ratna[!is.finite(ratna)] <- NA
  
  ut <- rowMeans(ratna, na.rm = TRUE)
  ut[k]
  ut_matrix[, replicate_num] <- ut
}

plot(design, ut_matrix[,1])
#lfire max utility is negative, -0.0014, optimal design =3


