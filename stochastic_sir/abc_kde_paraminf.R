rm(list=ls())
library(deSolve)
library(doParallel)
totalCores <- detectCores()
cluster <- makeCluster(totalCores[1]) 
registerDoParallel(cluster)
library(foreach)
library(truncnorm)
library(odin)

design<- seq(0,3,by=0.1)
num_replicates <- 1
ut_matrix <- matrix(NA, nrow = length(design), ncol = num_replicates)
S_init <- 490
I_init <- 10
R_init <- 0
dt <- 0.1
steps <- 30

start<-Sys.time()
for(replicate_num in 1:num_replicates){
  
  library(deSolve)
  library(doParallel)
  simulate_sir<- function(beta, gamma, S_init, I_init, R_init, dt, steps) {
    # Initialize vectors to store S, I, R
    S <- numeric(steps + 1)
    I <- numeric(steps + 1)
    R <- numeric(steps + 1)
    
    # Set initial values
    S[1] <- S_init
    I[1] <- I_init
    R[1] <- R_init
    for (i in 1:steps) {
      # Calculate increments with binomial distributions
      set.seed(123)
      delta_I <- rbinom(1, S[i], beta * I[i] / (S[i] + I[i] + R[i]))
      delta_R <- rbinom(1, I[i], gamma)
      
      # Update S, I, R using Euler's method
      S[i + 1] <- S[i] - delta_I
      I[i + 1] <- I[i] + delta_I - delta_R
      R[i + 1] <- R[i] + delta_R
    }
    
    
    # Return the result
    return(list(S = S, I = I, R = R))
  }
  
  
  
  
  cat("Replicate:", replicate_num, "\n")
  
  # Set a different seed for each replicate
  set.seed(replicate_num)
  
  # Priors
  beta = runif(1000, 0, 0.5) ##infection parameter
  gamma = runif(1000, 0, 0.5) ##recovery parameter
  parameters_values <- priorabc<-cbind(beta, gamma)
  
  z_recov_prior <- parameters_values
  sir_results <- lapply(1:nrow(priorabc), function(i) {
    as.data.frame(simulate_sir(priorabc[i, "beta"], priorabc[i, "gamma"], 
                               S_init, I_init, R_init, dt, steps))
  })
  design<- seq(0,3, by=0.1)
  
  trans <- vector("list", length = length(design))
  result_matrix <- matrix(NA, nrow = nrow(priorabc), ncol = 3)
  
  for (s in 1: length(design)) {
    for (i in 1:nrow(priorabc)){
      sir_res <- sir_results[[i]][s,]
      result_matrix[i, 1] <-sir_res[1,1]
      result_matrix[i, 2] <-sir_res[1,2]
      result_matrix[i, 3] <-sir_res[1,3]
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
      start<-Sys.time()
      cv.nn<-cv4abc(priorabc, trans[[k]][,1:3], nval=1000,
                     tols=c(0.09,0.1,0.2,0.3,0.4,0.5),
                     method="neuralnet",
             transf=c("logit", "logit"),
             logit.bounds=rbind(c(0,0.5), c(0,0.5)))
      stop<- Sys.time()
      timetaken<-stop-start
      summary(cv.nn)
      
      check<-abc(trans[[k]][i,1:3],priorabc, trans[[k]][,1:3],
                 tol=0.08,method="neuralnet", transf=c("logit", "logit"),
                 logit.bounds=rbind(c(0,0.5), c(0,0.5)))
      z.abc[[k]][[i]] <- check$adj.values[complete.cases(check$adj.values), ]
      
      set.seed(i)
      thetai <- matrix(NA, nrow = nrow(check$adj.values), ncol = 2)
      thetai[, 1] <- runif(nrow(check$adj.values), 0, 0.5)
      thetai[, 2] <- runif(nrow(check$adj.values), 0, 0.5)
      colnames(thetai) <- c("beta", "gamma")
      zi <- thetai[complete.cases(thetai), ]
      
      library(ks)
      set.seed(56)
      kde1<-kde(z.abc[[k]][[i]], density=TRUE)
      kde1p<- predict(kde1, x=z_recov_prior[i,])
      kde2<-kde(zi, density=TRUE)
      kde2p<- predict(kde2, x=z_recov_prior[i,])
      rati[k, i] <- kde1p/kde2p
    }
  }
  ratna <- log(rat)
  ratna[!is.finite(ratna)] <- NA
  
  ut <- rowMeans(ratna, na.rm = TRUE)
  ut_matrix[, replicate_num] <- ut
}
end<-Sys.time()
timetaken<-end-start
save(ut_matrix, file="utmatrixstochasticsirparaminf.RDa")

stopCluster(cluster)
library(ggplot2)
load("utmatrixstochasticsirparaminf.RDa")
mean_row <- rowMeans(ut_matrix, na.rm=TRUE)
sd_row <- apply(ut_matrix, 1, sd, na.rm=TRUE)
df <- data.frame(sequence = 1:nrow(ut_matrix), mean = mean_row, sd = sd_row)
pdf("sirparaminf.pdf", height=4, width=8)
ggplot(data=df, aes(x=design, y=mean)) + 
  labs(x = expression(paste(italic("d"))), 
       y = bquote(paste(italic(U[z](d)))), 
       color = "") +  theme_minimal() +
  geom_line() + theme(text=element_text(size=30))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha = 0.2)
dev.off()
design[which.max(df$mean)] #0.4
df$mean[which.max(df$mean)] #2.066
design[which.min(df$sd)]
df$sd[which.min(df$sd)]
df$sd[which.max(df$sd)]

data_df <- data.frame(Column1 = ut_matrix[, 1],
                      Column2 = ut_matrix[, 2],
                      Column3 = ut_matrix[, 3],
                      Column4 = ut_matrix[, 4])
library(RColorBrewer)
# Create a plot using ggplot2
library(ggplot2)
ggplot(data_df, aes(x = design)) +
  geom_line(aes(y = Column1, color = "Replicate 1"), size = 1) +
  geom_line(aes(y = Column2, color = "Replicate 2"), size = 1) +
  geom_line(aes(y = Column3, color = "Replicate 3"), size = 1) +
  geom_line(aes(y = Column4, color = "Replicate 4"), size = 1) +
  scale_color_manual(
    values = c(
      "Replicate 1" = "red",
      "Replicate 2" = "blue",
      "Replicate 3" = "green",
      "Replicate 4" = "purple"))+
  
  # "Replicate 6" = "pink"
  
  labs(
    x = expression(paste("Design, ", italic("d"))),
    y = bquote(paste("Utility, ", italic(U[z](d)))),
    color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(min(ut_matrix[, 1:4]), max(ut_matrix[, 1:4]) * 2.1)) +
  scale_y_continuous(
    breaks = seq(floor(min(ut_matrix[-1, 1:4])), ceiling(max(ut_matrix[-1, 1:4]) * 2.1), by = 0.002),
    labels = scales::number_format(accuracy = 0.001)
  )

design[which.max(ut)]
ut[which.max(ut)]

