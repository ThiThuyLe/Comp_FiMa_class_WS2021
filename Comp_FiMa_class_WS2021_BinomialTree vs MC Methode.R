# Jamie Lee Treber
# 568370
# Thema: Optionsbewertung: Binomialbaum Modell vs Monte Carlo-Methode

# Für das Erklärungsvideo bitte dem Link folgen:
# https://mediathek.htw-berlin.de/m/fa169e0553378d16cf0ec64d79045709339f0d98d9f8d63e96bf2637c96ead119ba73d7fec421611e3ba9698ee9d6a2989f7d137d5b3a5670a6995823fdee2c0


###################################################
###################################################

############## Binomialbaum Modell ################

###################################################
###################################################


# Zwei-Perioden-Binomialbaummodell

r <- 1/4 #Zinssatz, um aufzuzinsen
u <- 3/2 #up (Bernoulli - Zufallszahl)
d <- 1/2 #down (Bernoulli - Zufallszahl)
S0 <- 4  #Startwert zum Zeitpunkt 0
K <- 5   #Strike-Preis

S <- matrix(nrow = 3,ncol = 3,NA) #Matrix mit allen Preisen
S[1,1] = S0
S[1,2] = S[1,1]*u
S[2,2] = S[1,1]*d
S[1,3] = S[1,2]*u
S[2,3] = S[2,2]*u
S[3,3] = S[2,2]*d
S

#Werteprozess des Calls
V <- matrix(nrow = 3,ncol = 3,NA)
V[1,3] = max(S[1,3]-K,0)
V[2,3] = max(S[2,3]-K,0)
V[3,3] = max(S[3,3]-K,0)
V

# Risikoneutrales Wahrscheinlichkeitsmaß
q = (1+r-d)/(u-d)   #arbitragefreier Markt
q


V[1,2] = 1/(1+r)*( q*V[1,3] + (1-q)*V[2,3] ) 
V[2,2] = 1/(1+r)*( q*V[2,3] + (1-q)*V[3,3] )
V[1,1] = 1/(1+r)*( q*V[1,2] + (1-q)*V[2,2] )
V

#Werteprozess des Puts
V <- matrix(nrow = 3,ncol = 3,NA)
V[1,3] = max(K - S[1,3],0)
V[2,3] = max(K - S[2,3],0)
V[3,3] = max(K - S[3,3],0)
V


V[1,2] = 1/(1+r)*( q*V[1,3] + (1-q)*V[2,3] ) 
V[2,2] = 1/(1+r)*( q*V[2,3] + (1-q)*V[3,3] )
V[1,1] = 1/(1+r)*( q*V[1,2] + (1-q)*V[2,2] )
V



###############################################
###############################################

# n-Perioden Binomialbaummodell

start.time <- Sys.time()

# Hull, S. 476

r <- 0.1      #Zinssatz
T <- 5/12     #Laufzeit
N <- 5        #Aufteilung des Intervals in N Schritte
K <- 50       #Strike-Preis
dt <- T/N     #delta t - Laufzeitintervalle
S0 <- 50      #Preis des Underlying
sigma <- 0.4  #Volatilität

u <- exp(sigma*sqrt(dt)) #up
d <- 1/u                 #down
q <- ( exp(r*dt) - d  )/( u - d )   #Risikoneutrales W-Maß
dis <- exp(-r*dt)   #Diskontierung

# Matrizen deklarieren
S = matrix(nrow = N+1, ncol = N+1)
P = matrix(nrow = N+1, ncol = N+1)
C = matrix(nrow = N+1, ncol = N+1)
P.tmp = matrix(nrow = N+1, ncol = N+1)
C.tmp = matrix(nrow = N+1, ncol = N+1)
Exercise.P = matrix(nrow = N+1, ncol = N+1)
Exercise.C = matrix(nrow = N+1, ncol = N+1)

# Kurse für den Basiswert berechnen 
for (j in 1:(N+1)){
  for (i in (1:j)){
    S[i,j] <- S0*u^(j-i)*d^(i-1)
  }
}

# Auszahlungen bei Fälligkeit  
for (j in 1:(N+1)){
  P[j,N+1] <- max(K - S[j,N+1], 0)
  C[j,N+1] <- max(S[j,N+1] -K , 0)
  Exercise.P[j,N+1] <- ifelse(K - S[j,N+1] > 0,1,0)
  Exercise.C[j,N+1] <- ifelse(S[j,N+1] - K > 0,1,0)
}

# Rückwärtzberechnung der Optionswerte  

for (j in (N:1)){
  for (i in (1:j)){
    
    P.tmp[i,j] <- dis*(q*P[i,j+1] +(1-q)*P[i+1,j+1])   
    P[i,j] <- max(P.tmp[i,j], max(K - S[i,j],0) )  
    
    if(P[i,j] > P.tmp[i,j])
    { Exercise.P[i,j] <- 1 }
    else
    { Exercise.P[i,j] <- 0 }
    
    C.tmp[i,j] <- dis*(   q*C[i,j+1]  +(1-q)*C[i+1,j+1]   )   
    C[i,j] <- max(C.tmp[i,j], max(S[i,j] - K,0) )
    
    if( C[i,j] > C.tmp[i,j] ) 
    { Exercise.C[i,j] <- 1 }
    else
    { Exercise.C[i,j] <- 0 }
    
  }
}


if(N<6) {
  cat(" ", "\n")
  cat("Die Optionswerte P", "\n")
  print(round(P, digits = 2))
  P[1, 1]
} else P[1, 1]


if(N<6) {
  cat(" ", "\n")
  cat("Die Optionswerte C", "\n")
  print(round(C, digits = 2))
  C[1, 1]
} else C[1, 1]


######################################
end.time <- Sys.time()
time.takenNPBBM <- end.time - start.time
time.takenNPBBM


###################################################
###################################################

############## Monte Carlo-Methode ################

###################################################
###################################################

# https://www.r-bloggers.com/2020/12/pricing-of-european-options-with-monte-carlo/

start.time <- Sys.time()
# Call- und Put-Option für Monte-Carlo
call_put_mc<-function(nSim=1000000, dt, r, sigma, S0, K) {
  
  Z <- rnorm(nSim, mean=0, sd=1)
  WT <- sqrt(dt) * Z
  ST = S0*exp((r - 0.5*sigma^2)*dt + sigma*WT)
  
  # Preis und Standardabweichung für die Call-Option
  simulated_call_payoffs <- exp(-r*dt)*pmax(ST-K,0)
  price_call <- mean(simulated_call_payoffs)
  sterr_call <- sd(simulated_call_payoffs)/sqrt(nSim)
  
  # Preis und Standardabweichung für die Put-Option
  simulated_put_payoffs <- exp(-r*dt)*pmax(K-ST,0)
  price_put <- mean(simulated_put_payoffs)
  sterr_put <- sd(simulated_put_payoffs)/sqrt(nSim)
  
  
  output<-list(price_call=price_call, sterr_call=sterr_call, 
               price_put=price_put, sterr_put=sterr_put)
  return(output)
  
}

n <- 1000000   #Anzahl der Simulationen
dt <- 0.5      #delta t - Zeitintervalle
r <- 0.1       #Zinssatz
sigma <- 0.4   #Volatilität
S0 <- 50       #Startwert
K <- 50        #Strike-Preis


set.seed(1)
results<-call_put_mc(n, dt, r, sigma, S0, K)

results

###########################
end.time <- Sys.time()
time.takenMC <- end.time - start.time
time.takenMC


###################################################
###################################################

################### Vergleich #####################

###################################################
###################################################


###################################################
################ Ergebnisse #######################
###################################################


# N-Perioden Binomialbaum
resultsBBM <- c(C[1,1],P[1,1])
resultsBBM
# Call: 6.359546 
# Put: 4.488459

# Monte Carlo-Methode
results
# Call: 6.79244
# Put: 4.352614

###################################################
#################### Zeit #########################
###################################################

# N-Perioden Binomialbaum
time.takenNPBBM
#Time difference of 0.08200002 secs

# Monte Carlo-Methode
time.takenMC
#Time difference of 0.1932511 secs


