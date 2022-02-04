# ------------------------------------------------------------------------------
# Published in:     Comp_FiMa_class_WS2021
# ------------------------------------------------------------------------------
# Project title:    Trinomial tree computation 
#                   comparison of binomial and trinomial trees 
#                   calculating price of European call   
# ------------------------------------------------------------------------------
# Description:      calculating price of European Call with binomial and 
#                   trinomial tree, comparing convergence against Black-Scholes
#                   price
# ------------------------------------------------------------------------------
# Output:           Convergence of price of an European Call.png
# ------------------------------------------------------------------------------
# Author:           Patricia Ebert 
# ------------------------------------------------------------------------------
# Input:            S0: current stock price, X: strike price,   
#                   r: interest rate, T:  time to maturity,    
#                   sigma: volatility, N: number of steps 
# ------------------------------------------------------------------------------

rm(list = ls())

## functions for binomial and trinomial tree and Black-Scholes price 
binomialtree <- function(S0,X,r,T,sigma,N) {
  # compute constants
  dt <- T/N; v <- exp(-r * dt)
  u  <- exp(sigma * sqrt(dt)); d <- 1/u
  p  <- (exp(r * dt) - d) / (u - d)
  
  # initialize asset prices at maturity (period N)
  S <- numeric(N + 1)                                   
  for (j in 0:N) S[j+1] <- S0*u^j*d^(N-j)                
  
  # initialize option values at maturity (period N)
  C <- matrix(NA, nrow = N+1, ncol = N+1)
  for (j in 1:(N+1)) C[j,N+1] = pmax(S[j]-X, 0)
  
  # step back through the tree
  for (i in seq(N, 1, by = -1)){
    for(j in 1:i){
      C[j,i] <- v * ((1-p) * C[j,i+1] + p* C[j+1,i+1])
    }
  }
  C[1,1]
}
binomialtree(20, 20, 0.02, 5, 0.15, 1000)  

trinomialtree <- function(S0,X,r,T,sigma,N) {
  # compute constants
  dt <- T/N; v <- exp(-r * dt)
  u  <- exp(sigma * sqrt(2*dt)); d <- 1/u
  pu  <- ((exp((r*dt)/2) - exp(-sigma*sqrt(dt/2))) / (exp(sigma*sqrt(dt/2)) - exp(-sigma*sqrt(dt/2))))^2
  pd <- ((exp(sigma*sqrt(dt/2)) - exp((r*dt)/2)) / (exp(sigma*sqrt(dt/2)) - exp(-sigma*sqrt(dt/2))))^2
  
  # initialize asset prices at maturity (period N)
  S <- numeric(2*N + 1)                                      
  for (j in 0:(2*N)) S[j+1] <- S0*u^(-N+j)                   
  
  # initialize option values at maturity (period N)
  C <- matrix(NA, nrow = 2*N+1, ncol = N+1)
  for (j in 1:(2*N+1)) C[j,N+1] = pmax(S[j]-X, 0)
  
  # step back through the tree
  for (i in seq(N, 1, by = -1)){
    for(j in 1:(2*i-1)){
      C[j,i] <- v * (pd * C[j, i+1] + (1-pd-pu)*C[j+1, i+1] + pu*C[j+2,i+1])
    }
  }
  return (C[1,1])
}
trinomialtree(20, 20, 0.02, 5, 0.15, 1000)  

blackscholes <- function(S, X, rf, T, sigma) {
  d1 = (log(S/X) + (rf + sigma^2/2) * T)/(sigma * sqrt(T)) 
  d2 = d1 - sigma * sqrt(T)
  value = S * pnorm(d1) - X * exp(-rf * T) * pnorm(d2)
  value
}
blackscholes(20, 20, 0.02, 5, 0.15) 


## graphical representation of the convergence behaviour 

# price calculation for increasing N 
N <- seq(100, 1000, 100)
priceBinomial <- numeric(length(N))
priceTrinomial <- numeric(length(N))
for (i in 1:length(N)){
  steps <- N[i]
  priceBinomial[i]  <- binomialtree(20, 20, 0.02, 5, 0.15, steps)    
  priceTrinomial[i] <- trinomialtree(20, 20, 0.02, 5, 0.15, steps)
} 

# graphic 
plot(N, priceTrinomial, col = "blue", pch = 19,     
     main = "Convergence of price of an European Call \n for trinomial and binomial trees", 
     xlab = "number of steps", ylab = "price of European Call", 
     ylim = c(min(priceTrinomial, priceBinomial) - 0.001, max(priceTrinomial, priceBinomial) + 0.001))  
points(N, priceBinomial, col = "red", pch = 19)

bs <- blackscholes(20, 20, 0.02, 5, 0.15)
abline(h = bs, col = "darkgreen", lty = "dashed", lwd = 1.5)

legend("bottomright", 
       legend = c("Trinomial Tree", "Binomial Tree", "Black-Scholes"),
       lwd = 2, col = c("blue", "red", "darkgreen"), 
       lty = c(1,1,2))


## deviation from the BS-price
# absolute and relative error of trinomial and binomial tree  
relative_error_T <- round((abs(priceTrinomial - bs) / bs) * 100, 4) # in percent 
relative_error_B <- round((abs(priceBinomial - bs) / bs) * 100, 4) 
absolute_error_T <- round((abs(priceTrinomial - bs)), 4)
absolute_error_B <- round((abs(priceBinomial - bs)), 4)

df <- data.frame(N, bs, priceTrinomial, priceBinomial, absolute_error_T, absolute_error_B, relative_error_T, relative_error_B)
colnames(df) <- c("steps", "Black-Scholes", "trinomial", "binomial", "abs. error binomial", "abs. error trinomial", "rel. error trinomial", "rel. error binomial")
