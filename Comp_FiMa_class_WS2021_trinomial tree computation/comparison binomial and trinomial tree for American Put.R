# ------------------------------------------------------------------------------
# Published in:     Comp_FiMa_class_WS2021
# ------------------------------------------------------------------------------
# Project title:    Trinomial tree computation 
#                   comparison of binomial and trinomial trees 
#                   calculating price of American Put    
# ------------------------------------------------------------------------------
# Description:      calculating price of American Put with binomial and 
#                   trinomial tree, comparing convergence 
# ------------------------------------------------------------------------------
# Output:           Convergence of price of an American Put.png
# ------------------------------------------------------------------------------
# Author:           Patricia Ebert 
# ------------------------------------------------------------------------------
# Input:            S0: current stock price, X: strike price,   
#                   r: interest rate, T:  time to maturity,    
#                   sigma: volatility, N: number of steps 
# ------------------------------------------------------------------------------

rm(list = ls())

## functions for binomial and trinomial tree 
binomialtreeAP <- function(S0,K,r,T,sigma,N){
  # compute constants
  dt <- T/N; dis <- exp(-r*dt)
  u <- exp(sigma * sqrt(dt)); d <- 1/u 
  q <- (exp(r*dt)-d)/(u-d)
  
  # value of underlying  
  S = matrix(nrow = N+1, ncol = N+1)
  for (j in 1:(N+1)){
    for (i in (1:j)){
      S[i,j] <- S0*u^(j-i)*d^(i-1)
    }
  }
  
  # option values at maturity   
  P = matrix(nrow = N+1, ncol = N+1)
  Exercise.P = matrix(nrow = N+1, ncol = N+1)
  for (j in 1:(N+1)){
    P[j,N+1] <- max(K - S[j,N+1], 0)
    Exercise.P[j,N+1] <- ifelse(K - S[j,N+1] > 0,1,0)   
    # exercise = 1 if option value > 0  
  }
  
  # step back through the tree 
  P.tmp = matrix(nrow = N+1, ncol = N+1)
  for (j in (N:1)){
    for (i in (1:j)){
      
      P.tmp[i,j] <- dis*(q*P[i,j+1]+(1-q)*P[i+1,j+1])   
      P[i,j] <- max(P.tmp[i,j], max(K - S[i,j],0))  
      
      if(P[i,j] > P.tmp[i,j])
      { Exercise.P[i,j] <- 1 }     
      else
      { Exercise.P[i,j] <- 0 }
    }
  }
  P[1,1]
}
binomialtreeAP(50, 50, 0.1, 5/12, 0.4, 1000) 


trinomialtreeAP <- function(S0,K,r,T,sigma,N){
  # compute constants
  dt <- T/N; dis <- exp(-r*dt)
  u  <- exp(sigma * sqrt(2*dt)); d <- 1/u
  pu  <- ((exp((r*dt)/2) - exp(-sigma*sqrt(dt/2))) / (exp(sigma*sqrt(dt/2)) - exp(-sigma*sqrt(dt/2))))^2
  pd <- ((exp(sigma*sqrt(dt/2)) - exp((r*dt)/2)) / (exp(sigma*sqrt(dt/2)) - exp(-sigma*sqrt(dt/2))))^2
  
  
  # value of underlying  
  S = matrix(nrow = 2*N+1, ncol = N+1)
  for (j in 1:(N+1)){
    for (i in (1:(2*j-1))){
      S[i,j] <- S0*u^(j-i)
    }
  }
  
  # option values at maturity   
  P = matrix(nrow = 2*N+1, ncol = N+1)
  Exercise.P = matrix(nrow = 2*N+1, ncol = N+1)
  for (j in 1:(2*N+1)){
    P[j,N+1] <- max(K - S[j,N+1], 0)
    Exercise.P[j,N+1] <- ifelse(K - S[j,N+1] > 0,1,0)   
    # exercise = 1 if option value > 0  
  }

  # step back through the tree 
  P.tmp = matrix(nrow = 2*N+1, ncol = N+1)
  for (j in (N:1)){
    for (i in (1:(2*j-1))){
      
      P.tmp[i,j] <- dis*(pu*P[i,j+1] + (1-pu-pd)*P[i+1,j+1] + pd*P[i+2, j+1])   
      P[i,j] <- max(P.tmp[i,j], max(K - S[i,j],0))  
      
      if(P[i,j] > P.tmp[i,j])
      { Exercise.P[i,j] <- 1 }     
      else
      { Exercise.P[i,j] <- 0 }
    }
  }
  P[1,1]
}
trinomialtreeAP(50, 50, 0.1, 5/12, 0.4, 1000) 


## graphical representation of the convergence behaviour 
# price calculation for increasing N 
N <- seq(100, 1000, 100)
priceBinomial <- numeric(length(N))
priceTrinomial <- numeric(length(N))
for (i in 1:length(N)){
  steps <- N[i]
  priceBinomial[i]  <- binomialtreeAP(50, 50, 0.1, 5/12, 0.4, steps)
  priceTrinomial[i] <- trinomialtreeAP(50, 50, 0.1, 5/12, 0.4, steps)
} 

# graphic 
plot(N, priceTrinomial, col = "blue", pch = 19, cex = 0.9,     
     main = "Convergence of price of an American Put \n for trinomial and binomial trees", 
     xlab = "number of steps", ylab = "price of American Put", 
     ylim = c(min(priceTrinomial, priceBinomial) - 0.0001, max(priceTrinomial, priceBinomial) + 0.001))  
points(N, priceBinomial, col = "red", pch = 19, cex = 0.9)

legend("bottomright", 
       legend = c("Trinomial Tree", "Binomial Tree"),
       lwd = 2, col = c("blue", "red"), 
       lty = c(1,1))


## difference between price of trinomial and binomial tree 
difference <- round((abs(priceTrinomial - priceBinomial)), 4)

df <- data.frame(N, priceTrinomial, priceBinomial, difference)
colnames(df) <- c("steps","trinomial", "binomial", "deviation")

