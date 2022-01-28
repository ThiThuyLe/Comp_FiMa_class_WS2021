################################################
# comparison of binomial and trinominal trees 
# calculating price of European call   
################################################

# S0: current stock price, X: strike price,   r: interest rate
# T:  time to maturity,    sigma: volatility, N: number of steps   

rm(list = ls())

S0 = 20; X = 20; r = 0.02 
T = 5; sigma = 0.15; N = 2 

Binomialbaum <- function(S0,X,r,T,sigma,N) {
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
Binomialbaum(20, 20, 0.02, 5, 0.15, 1000)  

Trinomialbaum <- function(S0,X,r,T,sigma,N) {
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
Trinomialbaum(20, 20, 0.02, 5, 0.15, 1000)  

blackscholes = function(S, X, rf, T, sigma) {
  d1 = (log(S/X) + (rf + sigma^2/2) * T)/(sigma * sqrt(T)) 
  d2 = d1 - sigma * sqrt(T)
  value = S * pnorm(d1) - X * exp(-rf * T) * pnorm(d2)
  value
}
blackscholes(20, 20, 0.02, 5, 0.15) # 3.59934 


## graphical representation of the convergence behaviour 

# price calculation for increasing N 
N <- seq(100, 1000, 100)
priceBinomial <- numeric(length(N))
priceTrinomial <- numeric(length(N))
for (i in 1:length(N)){
  steps <- N[i]
  priceBinomial[i]  <- Binomialbaum(20, 20, 0.02, 5, 0.15, steps)
  priceTrinomial[i] <- Trinomialbaum(20, 20, 0.02, 5, 0.15, steps)
} 

plot(N, priceTrinomial, col = "blue", pch = 19,     
     main = "Convergence of Price of an European Call \n for Trinomial and Binomial Trees", 
     xlab = "number of steps", ylab = "Price of European Call", 
     ylim = c(min(priceTrinomial, priceBinomial) - 0.001, max(priceTrinomial, priceBinomial) + 0.001))  
points(N, priceBinomial, col = "red", pch = 19)

bs <- blackscholes(20, 20, 0.02, 5, 0.15)
abline(h = bs, col = "darkgreen", lty = "dashed", lwd = 1.5)

legend("bottomright", 
       legend = c("Trinomial Tree", "Binomial Tree", "Black Scholes"),
       lwd = 2, col = c("blue", "red", "darkgreen"), 
       lty = c(1,1,2))
