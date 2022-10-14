# HEADER ####
# Poisson-based sampling transect length
# PD
# 20 Oct 2021, updated 10 and 23 Nov 2021

# The script outputs a table to justify how long a sampling transect should be 
# for a given expected (average) number of plants.
# The script does not make an optimization - the range of transect lengths has to be tried out.
# Assumption: the distribution of plants is such that samples have Poisson distribution.
# Settings: 
# - expected (average) number of plants (seeds etc.) per unit length of row
# - transect length
# - number of transects
# - number of simulations - optionally
# Output: a csv table of p values - how likely the count is to depart from global average by no more than a predefined tolerance. 

# SETTINGS ----------------------------------------------------------------------------------------------------
mu <- 6 # Global means (expected number of plants per unit length, e.g. from drill settings or from field survey)
# (56 - New Hopyard 2018; 10.3 - New Hopyard 2020)
N <- 12*8 # Total number of transects 
#L <- seq(from=0.10, by=0.05, length.out=110) # vector of transect lengths
L <- seq(from=1, by=0.10, length.out=110) # vector of transect lengths

d <- seq(from=0.01, by=0.01, length.out=10) # vector of tolerances for the mean, one-sided ("departure from the means by +/-d")
# (16*22 - New Hopyard one treatment)
nsim <- 10000 # Number of simulations - adjust optionally

# RUN ----------------------------------------------------------------------------------------------------------
mt <- mu*L # expected number of plants per transect
m <- matrix(NA, nrow=length(L), ncol=nsim) # dummy matrix for simulated means
p.value <- matrix(NA, nrow=length(L), ncol=length(d)) # dummy vector for p.value

# Simulations for a range of transect lengths
for(j in 1:length(L)){ # for every transect length in the range
  r <- rpois(n = N*nsim, lambda = mt[j])   # simulation nsim times, with N transects every time
  r <- matrix(r, ncol=nsim) # ... folded into a matrix where every column stands for a separate simulation
  m[j,] <- colMeans(r) # calculating the mean for every simulation
}

# Calculate p-value 
for(j in 1:length(L)){ # for every transect length
  for(k in 1:length(d)){ # for every tolerance
    # find the simulations where the mean departed from the expected value by more than the tolerance
    # count them and relate to the overall number of simulations (calculate probability)
#   p.value[j,k] <- sum(m[j,] > (1-d[k])*mt[j] & m[j,] < (1+d[k])*mt[j]) / nsim    # optionally for 1-p.value
    p.value[j,k] <- sum(m[j,] < (1-d[k])*mt[j] | m[j,] > (1+d[k])*mt[j]) / nsim    # optionally for p.value
  }
}
p.value <- round(p.value,4)

# Save result as csv for further reference -----------------------------------------------------------------------
result <- data.frame(mu, N, L, p.value)
colnames(result) <- c("Unit density", "Number of transects", "Transect length", paste0("_", d*100, "%"))
fld.path <- 'output'
filename <- paste0(fld.path, "/Poisson transect length", " N-", N, " mu-", mu, ".csv")
write.csv(result, filename, row.names = F)
#rm(list=setdiff(ls(), c("result")))
