rm(list=ls())
library(deSolve)
library(doParallel)
totalCores = detectCores()
cluster <- makeCluster(totalCores[1]) #, outfile="" 
registerDoParallel(cluster)
library(foreach)
library(truncnorm)

# Define the FitzHugh-Nagumo model and compute spikes
fitzhugh_nagumo<- function(time, y0, parms) {
  n <- length(time) - 1
  dt <- diff(time)
  y <- matrix(NA, nrow = 2, ncol = n + 1)
  y[, 1] <- y0
  
  with(as.list(parms), {
    for (i in 1:n) {
      # Generate Brownian motion increments
      sigma <- 0.1
      set.seed(123)
      dB <- rnorm(1, mean = 0, sd = sqrt(dt[i]) * sigma)
      
      # Stochastic differential equations with Brownian motion
      du_dt <- parms[1] * (y[1, i] - y[1, i]^3 / 3 + y[2, i] + parms[4]) 
      dv_dt <- -(1 / parms[4]) * (y[1, i] - parms[2] + parms[3] * y[2, i]) + dB
      
      # Print diagnostic information
      #cat("Step:", i, "\n")
      #cat("Time:", time[i], "\n")
      #cat("du_dt:", du_dt, "\n")
      #cat("dv_dt:", dv_dt, "\n")
      
      # Euler-Maruyama integration
      y[, i + 1] <- y[, i] + c(du_dt, dv_dt) * dt[i]
      #cat("y:", y[, i + 1], "\n\n")
    }
    return(t(y))
  })
  
}
compute_spikes <- function(parms, uspike, time, y_initial) {
  fitzhugh_nagumo<- function(time, y0, parms) {
    n <- length(time) - 1
    dt <- diff(time)
    y <- matrix(NA, nrow = 2, ncol = n + 1)
    y[, 1] <- y0
    
    with(as.list(parms), {
      for (i in 1:n) {
        # Generate Brownian motion increments
        sigma <- 0.1
        set.seed(123)
        dB <- rnorm(1, mean = 0, sd = sqrt(dt[i]) * sigma)
        
        # Stochastic differential equations with Brownian motion
        du_dt <- parms[1] * (y[1, i] - y[1, i]^3 / 3 + y[2, i] + parms[4]) 
        dv_dt <- -(1 / parms[1]) * (y[1, i] - parms[2] + parms[3] * y[2, i]) + dB
        
        # Print diagnostic information
        #cat("Step:", i, "\n")
        #cat("Time:", time[i], "\n")
        #cat("du_dt:", du_dt, "\n")
        #cat("dv_dt:", dv_dt, "\n")
        
        # Euler-Maruyama integration
        y[, i + 1] <- y[, i] + c(du_dt, dv_dt) * dt[i]
        #cat("y:", y[, i + 1], "\n\n")
      }
      return(t(y))
    })
    
  }
  
  set.seed(67)
  solution <- fitzhugh_nagumo(time, y_initial, parms)
  u_t=solution[,2]
  # Compute spike rate and average spike duration
  epsilon <- 1e-6 # Set the interval for looking back in time (epsilon)
  # Initialization
  t_spike <- NULL
  start_spikes <- NULL
  end_spikes <- NULL
  
  # Variable to keep track if we are currently in a spike
  in_spike <- FALSE
  
  # Detecting spikes
  for (t in 2:length(u_t)) {
    if (!in_spike) {
      if (uspike <= u_t[t] && uspike > u_t[max(1, t - epsilon)]) {
        t_spike <- c(t_spike, t) # append the spike time
        start_spikes <- c(start_spikes, t)
        in_spike <- TRUE
      }
    } else {
      if (u_t[t] < uspike) {
        end_spikes <- c(end_spikes, t - 1) # append the end time of the spike
        in_spike <- FALSE
      }
    }
  }
  
  if (in_spike) {
    end_spikes <- c(end_spikes, length(u_t))
  }
  
  num_spikes <- length(t_spike)
  
  spike_duration <- mean(end_spikes - start_spikes)
  
  return(list(spike_rate = num_spikes/0.2, spike_duration = spike_duration))
}

# Set parameters and initial conditions
uspike <- 0.5
time <- seq(0.1, 200, length=1000)
y_initial <- c(u = 0, v = 0)

#design variables: zeta
zetarange<- seq(from=0, to=0.8, length=30)

#generate 10k samples for ABC and data
set.seed(123)
priorabc<-cbind(rtruncnorm(1000, a=0,b=1, mean=0.4, sd=0.3),
                rtruncnorm(1000, a=0,b=1, mean=0.4, sd=0.4))
abcresults=list()
for(k in 1:length(zetarange)){
  abcresults[[k]]=matrix(NA, nrow=1000, ncol=4)
  for(i in 1:nrow(priorabc)){
    parms=c(gamma = 3.0, theta0 = priorabc[i,1],
            theta1 = priorabc[i,2], zeta = zetarange[k])
    cs=compute_spikes(parms, uspike, time, y_initial)
    if (!is.na(cs$spike_rate) && !is.na(cs$spike_duration) && cs$spike_duration != 0) {
      abcresults[[k]][i, 1] <- priorabc[i,1]
      abcresults[[k]][i, 2] <- priorabc[i,2]
      abcresults[[k]][i, 3] <- cs$spike_rate
      abcresults[[k]][i, 4] <- cs$spike_duration
    }
  }
}

num_replicates <- 5
ut_matrix <- matrix(NA, nrow = length(zetarange), ncol = num_replicates)


for(replicate_num in 1:num_replicates){
  
  cat(" Starting Replicate:", replicate_num, " at ", Sys.time(), "\n")   
  # Set a different seed for each replicate
  set.seed(replicate_num)
  theta0=rtruncnorm(1000, a=0,b=1, mean=0.4, sd=0.3)
  theta1=rtruncnorm(1000, a=0,b=1, mean=0.4, sd=0.4)
  parameters_values<-cbind(theta0,theta1)
  
  #goal: parameter inference
  #z_prior<- parameters_values
  cur<-0.2
  uspike <- 0.5
  z_prior<-c()
  for(l in 1:nrow(parameters_values)){
    parms=c(gamma = 3.0, theta0 = parameters_values[l,1],
            theta1 = parameters_values[l,2], zeta =cur)
    solution <- compute_spikes(parms, uspike, time, y_initial)
    z_prior[l]<- solution$spike_rate
  }
  
  spike_results=list()
  for(k in 1:length(zetarange)){
    spike_results[[k]]=matrix(NA, nrow=1000, ncol=4)
    for(i in 1:nrow(parameters_values)){
      parms=c(gamma = 3.0, theta0 = theta0[i],
              theta1 = theta1[i], zeta = zetarange[k])
      cs=compute_spikes(parms, uspike, time, y_initial)
      if (!is.na(cs$spike_rate) && !is.na(cs$spike_duration) && cs$spike_duration != 0) {
        spike_results[[k]][i, 1] <- theta0[i]
        spike_results[[k]][i, 2] <- theta1[i]
        spike_results[[k]][i, 3] <- cs$spike_rate
        spike_results[[k]][i, 4] <- cs$spike_duration
      }
    }
  }
  
  rat_list <- vector("list", length(zetarange))
  
  for(k in 2:length(zetarange)){
    cat(" Starting Replicate:", replicate_num, "with deisgn point",k,
        "at", print(Sys.time()), "\n")
    cc <- spike_results[[k]][complete.cases(spike_results[[k]]), ]
    # ccabc<-abcresults[[k]][complete.cases(abcresults[[k]]), ]
    rat_list[[k]] <- foreach(i = 1:nrow(cc), .combine = cbind) %dopar% {
      #for(i in 1:nrow(cc)){
      library(deSolve)
      library(truncnorm)
      library(abc)
      check<- abc(cc[i,3:4],cc[,1:2],cc[,3:4],method="rejection",
                  tol=0.05) #,
      #method="neuralnet", transf=c("logit", "logit"),
      #        logit.bounds=matrix(c(0,0,1,1), nrow=2, ncol=2))
      
      check$adj.values<-check$unadj.values
      check$adj.values <- check$adj.values[complete.cases(check$adj.values), ]
      
      
      if (nrow(check$adj.values) != 0) {
        
        fitzhugh_nagumo<- function(time, y0, parms) {
          n <- length(time) - 1
          dt <- diff(time)
          y <- matrix(NA, nrow = 2, ncol = n + 1)
          y[, 1] <- y0
          
          with(as.list(parms), {
            for (i in 1:n) {
              # Generate Brownian motion increments
              sigma <- 0.1
              set.seed(123)
              dB <- rnorm(1, mean = 0, sd = sqrt(dt[i]) * sigma)
              
              # Stochastic differential equations with Brownian motion
              du_dt <- parms[1] * (y[1, i] - y[1, i]^3 / 3 + y[2, i] + parms[4]) 
              dv_dt <- -(1 / parms[1]) * (y[1, i] - parms[2]+ parms[3] * y[2, i]) + dB
              
              # Print diagnostic information
              #cat("Step:", i, "\n")
              #cat("Time:", time[i], "\n")
              #cat("du_dt:", du_dt, "\n")
              #cat("dv_dt:", dv_dt, "\n")
              
              # Euler-Maruyama integration
              y[, i + 1] <- y[, i] + c(du_dt, dv_dt) * dt[i]
              #cat("y:", y[, i + 1], "\n\n")
            }
            return(t(y))
          })
          
        }
        cur<-0.2
        uspike<-0.5
        z.abc<-c()
        for(l in 1:nrow(check$adj.values)){
          parms=c(gamma = 3.0, theta0 = check$adj.values[l,1],
                  theta1 = check$adj.values[l,2], zeta =cur)
          solution <- compute_spikes(parms, uspike, time, y_initial)
          z.abc[l]<- solution$spike_rate
        }
        
        set.seed(i)
        theta0i=rtruncnorm(nrow(check$adj.values), a=0,b=1, mean=0.4, sd=0.3)
        theta1i=rtruncnorm(nrow(check$adj.values), a=0,b=1, mean=0.4, sd=0.4)
        thetai<-cbind(theta0i, theta1i)
        zi<-c()
        cur<-0.2
        uspike<-0.5
        for(l in 1:nrow(check$adj.values)){
          parms=c(gamma = 3.0, theta0 = theta0i[l],
                  theta1 = theta1i[l], zeta =cur)
          solution <- compute_spikes(parms, uspike, time, y_initial)
          zi[l]<- solution$spike_rate
        }
        library(densratio)
        set.seed(89)
        densratio_obj <- densratio::densratio(z.abc, zi,method="RuLSIF", alpha=0.5,
                                              sigma = 10 ^ seq(-5, 10, length.out = 50),
                                              lambda = 10 ^ seq(-5, 10, length.out = 50))
        pp<-t(as.matrix(cc[i,1:2]))
        densratio_obj$compute_density_ratio(z_prior)
      }
    }
  }
  ut<-c()
  for(k in 2:length(zetarange)){
    df<-rat_list[[k]]
    df_nonzero <- df[df>1]  # Remove elements equal to 0
    ut[k] <- mean(log(df_nonzero), na.rm=TRUE)
  }
  ut_matrix[, replicate_num] <- ut
}
save(ut_matrix, file="fhnreplicatesgoal.RDa") # for uspike 0.5
library(ggplot2)
design<- seq(from=0, to=0.8, length=30)
mean_row <- rowMeans(ut_matrix, na.rm=TRUE)
sd_row <- apply(ut_matrix, 1, sd)
df <- data.frame(sequence = 1:nrow(ut_matrix), mean = mean_row, sd = sd_row)
pdf("fhngoal.pdf", height=4, width=8)
ggplot(data=df, aes(x=zetarange, y=mean)) + 
  labs(x = expression(paste(italic("d"))), 
       y = bquote(paste(italic(U[z](d)))), 
       color = "") +  theme_minimal() +
  geom_line() + theme(text=element_text(size=30))+
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), alpha = 0.2)
dev.off()
zetarange[which.max(df$mean)] #0.772
df$mean[which.max(df$mean)] # 0.768
zetarange[which.min(df$sd)]
df$sd[which.min(df$sd)]
df$sd[which.max(df$sd)]

# Create a data frame
data_df <- data.frame(Column1 = ut_matrix[, 1],
                      Column2 = ut_matrix[, 2],
                      Column3 = ut_matrix[, 3],
                      Column4 = ut_matrix[, 4],
                      Column5 = ut_matrix[, 5])

# Create a plot using ggplot2
library(ggplot2)
ggplot(data_df, aes(x = zetarange)) +
  geom_line(aes(y = Column1, color = "Replicate 1"), size = 1) +
  geom_line(aes(y = Column2, color = "Replicate 2"), size = 1) +
  geom_line(aes(y = Column3, color = "Replicate 3"), size = 1) +
  geom_line(aes(y = Column4, color = "Replicate 4"), size = 1) +
  geom_line(aes(y = Column5, color = "Replicate 5"), size = 1) +
  scale_color_manual(values = c("Replicate 1" = "red",
                                "Replicate 2" = "blue",
                                "Replicate 3" = "green",
                                "Replicate 4" = "purple",
                                "Replicate 5" = "orange")) +
  labs(x = expression(paste("Design, ", italic("d"))), 
       y = bquote(paste("Utility, ", italic(U[z](d)))), 
       color = "") +  theme_minimal() +
  theme(legend.position = "top") +
  coord_cartesian(ylim = c(min(ut_matrix[-1, 1:5]), max(ut_matrix[-1, 1:5]) ))
