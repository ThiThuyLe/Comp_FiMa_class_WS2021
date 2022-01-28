#####################################
# Trinomial processes 
#####################################
# https://github.com/QuantLet/SFE/blob/master/QID-3167-SFEtrinomp/README.md

rm(list = ls(all = TRUE))

# declaration of variables  
pu = 0.4       # probability of up-movement  
pd = 0.25      # probability of down-movement  
n  = 100       # number of steps 
k  = 6         # number of paths 
u  = 1         # up-movement 
d  = 1.1       # down-movement 
t  = 1:n       

# pseudo random numbers (uniform)
set.seed(2022)
z = matrix(runif(n*k), n, k)
z = (-d) * (z < pd) + u*(z > (1-pu))
head(z)  
  # for z < pd z = -d 
  # for z > 1-pu z = u 
  # for pd < z < 1-pu = 0 
x = apply(z, MARGIN = 2, FUN = cumsum)
  # accumulate z column-wise
head(x)


# trend und upper und lower confidence band 
trend = t * (pu * u - pd * d)
std   = sqrt(t * (pu * (1-pu)*u^2 + pd * (1-pd) * d^2 + 2*pu*pd*u*d)) 
s1    = trend + 2*std    # upper confidence band 
s2    = trend - 2*std    # lower confidence band 

# plot of trinomial process 
matplot(x, main = paste("Trinomial process with pu =", pu, ", pd =", pd, ", pm =", 1 - pu - pd, 
        "\n and u =", u, ", d =", round(d,2)),
        xlab = "Time", ylab = "Process", type = "l", lwd = 2.5, 
        lty = 1, ylim = c(min(x, s2), max(x, s1)), col = 2:(k + 1))
lines(trend)
lines(s1, lty = 2)
lines(s2, lty = 2)




