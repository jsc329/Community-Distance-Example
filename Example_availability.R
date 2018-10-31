

model{
  
  beta0 ~ dunif(-10, 10)
  
  #covs on availability
  logit(p.a) <- beta0
  
  # For time period 1
  pi.pa[1] <- 1 - pow(p.a, int_len[1])
  pi.pa.c[1] <- pi.pa[1] / phi
  
  for (k in 2:K){ # Cycle through the remaining time periods
    # Essentially the first "pow" sums up all the time periods
    # up to one before the period of interest, since we want to know what the
    # chance is of NOT detecting an individual in the periods prior
    # multiplied by the chance of detecting a bird in the current period k
    pi.pa[k] <- pow(p.a, sum(int_len[1:(k-1)])) * (1 - pow(p.a, int_len[k]))
    pi.pa.c[k] <- pi.pa[k] / phi
  }
  
  phi <- sum(pi.pa) # Probability ever available at a site, based on singing rate
  
  for(j in 1:nobs){
    
    tint[j] ~ dcat(pi.pa.c[]) # Spot with actual data
    
  }
  
  
}