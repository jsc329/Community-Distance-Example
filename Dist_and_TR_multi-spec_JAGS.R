
model {

# This is taken from pg. 476 in Kery and Royle 2015
  
  #Hypers
  mu_b.a0 ~ dnorm(0, 0.01)
  mu_b.a1 ~ dnorm(0, 0.01)
  mu_b.a2 ~ dnorm(0, 0.01)
  mu_a0 ~ dnorm(0, 0.01)
  mu_a1 ~ dnorm(0, 0.01)
  mu_a2 ~ dnorm(0, 0.01)
  mu_b0 ~ dnorm(0, 0.01)
  mu_b1 ~ dnorm(0, 0.01)
  #mu_b2 ~ dnorm(0, 0.01)
  
  tau_b.a0 ~ dgamma(0.1, 0.1)
  tau_b.a1 ~ dgamma(0.1, 0.1)
  tau_b.a2 ~ dgamma(0.1, 0.1)
  tau_a0 ~ dgamma(0.1, 0.1)
  tau_a1 ~ dgamma(0.1, 0.1)
  tau_a2 ~ dgamma(0.1, 0.1)
  tau_b0 ~ dgamma(0.1, 0.1)
  tau_b1 ~ dgamma(0.1, 0.1)
  #tau_b2 ~ dgamma(0.1, 0.1)
  

  for (i in 1:nspec){
    # Priors for parameters
    beta.a0[i] ~ dnorm(mu_b.a0, tau_b.a0)
    alpha0[i] ~ dnorm(mu_a0, tau_a0)
    alpha1[i] ~ dnorm(mu_a1, tau_a1)
    alpha2[i] ~ dnorm(mu_a2, tau_a2)
    
    beta.a1[i] ~ dnorm(mu_b.a1, tau_b.a1)
    beta.a2[i] ~ dnorm(mu_b.a2, tau_b.a2)
    beta0[i] ~ dnorm(mu_b0, tau_b0)
    beta1[i] ~ dnorm(mu_b1, tau_b1)
    #beta2[i] ~ dnorm(mu_b2, tau_b2)

    
for (s in 1:nsites){
	for (t in 1:nvisits){

		#covariates on perceptibility
		log(sigma[s,t,i]) <- alpha0[i] + alpha1[i]*time[s,t] + alpha2[i]*day[s,t]
		#covs on availability
		logit(p.a[s,t,i]) <- beta.a0[i] + beta.a1[i]*time[s,t] + beta.a2[i]*day[s,t]
	

		# Distance sampling detection model
		for (b in 1:nD){
			log(g[b,s,t,i]) <- -mdpts[b]*mdpts[b] / (2*sigma[s,t,i]*sigma[s,t,i])
			f[b,s,t,i] <- (2*mdpts[b]*delta) / (B*B)
			pi.pd[b,s,t,i] <- g[b,s,t,i]*f[b,s,t,i]
			pi.pd.c[b,s,t,i] <- pi.pd[b,s,t,i] / pdet[s,t,i]
		}
		
	  # pdet is detection probability based on perception (AKA distance to bird)
		pdet[s,t,i] <- sum(pi.pd[,s,t,i])
	
		# Time-removal detection probabilities
		
		# For time period 1
		  pi.pa[1,s,t,i] <- 1 - pow(p.a[s,t,i], int_len[1])
		  pi.pa.c[1,s,t,i] <- pi.pa[1,s,t,i] / phi[s,t,i]
		  
		for (k in 2:K){ # Cycle through the remaining time periods
		  # Essentially the first "pow" sums up all the time periods
		  # up to one before the period of interest, since we want to know what the
		  # chance is of NOT detecting an individual in the periods prior
		  # multiplied by the chance of detecting a bird in the current period k
		  pi.pa[k,s,t,i] <- pow(p.a[s,t,i], sum(int_len[1:(k-1)])) * (1 - pow(p.a[s,t,i], int_len[k]))
			pi.pa.c[k,s,t,i] <- pi.pa[k,s,t,i] / phi[s,t,i]
		}

		phi[s,t,i] <- sum(pi.pa[,s,t,i]) # Probability ever available at a site, based on singing rate

    #Abundance model

		# marginal detection probability is the product of perceptibility (pdet) from distance data
		# and singing rate (phi) AKA. whether a bird is singing or visible during a survey
		# There's another peice too, movement during a survey, but we're not going
		# to deal with that today
		pmarg[s,t,i] <- pdet[s,t,i] * phi[s,t,i]
		
		# Spot with actual data (n is the number of birds of species i, during survey t, at site s
		n[s,t,i] ~ dbin(pmarg[s,t,i], M[s,i]) 
	
  } # End of the visit level loop

  # M is a latent variable that represents the "true" number of birds of species i, at site s
  # during the entire primary period (across all surveys)
	M[s,i] ~ dpois(lambda[s,i])
	# Here's where we model abundance relative to whether a point is in an OHV use area
	log(lambda[s,i]) <- beta0[i] + beta1[i]*ohvind[s]
	
	# GOF for abundance
	Mnew[s,i] ~ dpois(lambda[s,i])
	
	# Residuals
	Mresid[s,i] <- pow(sqrt(M[s,i]) - sqrt(lambda[s,i]), 2)
	Mnewresid[s,i] <- pow(sqrt(Mnew[s,i]) - sqrt(lambda[s,i]), 2)
	
} # End of the site level loop
    
  # Sum over all sites
  Mres_spec[i] <- sum(Mresid[1:nsites, i])
  Mnewres_spec[i] <- sum(Mnewresid[1:nsites, i])
	
	# Derived things
	Mtot[i] <- sum(M[,i])
	MtotOHV[i] <- sum(M[OHVsites, i])
	MtotnoOHV[i] <- sum(M[noOHVsites, i])

	PDETmean[i] <- mean(pdet[,,i])
	PHImean[i] <- mean(phi[,,i])

}#End of the species level loop
  
  # Sum over all species to get an omnibus statistic
  Mres_all <- sum(Mres_spec[1:nspec])
  Mnewres_all <- sum(Mnewres_spec[1:nspec])
    
# Conditional model for categorical covariates
  for(j in 1:nobs){
      dclass[j] ~ dcat(pi.pd.c[,site[j], visit[j], species[j] ]) # Spot with actual data
      tint[j] ~ dcat(pi.pa.c[,site[j], visit[j], species[j] ]) # Spot with actual data
      
      # GOF Section
      dclassnew[j] ~ dcat(pi.pd.c[,site[j], visit[j], species[j] ])
      tintnew[j] ~ dcat(pi.pa.c[,site[j], visit[j], species[j] ])
      
      dist_resid[j] <- pow(1- sqrt(pi.pd.c[dclass[j], site[j], visit[j], species[j] ]), 2)
      dist_resid_new[j] <- pow(1- sqrt(pi.pd.c[dclassnew[j], site[j], visit[j], species[j] ]), 2)
      
      tint_resid[j] <- pow(1- sqrt(pi.pa.c[tint[j], site[j], visit[j], species[j] ]), 2)
      tint_resid_new[j] <- pow(1- sqrt(pi.pa.c[tintnew[j], site[j], visit[j], species[j] ]), 2)
  }
  
  # GOF section for perceptability and availability
  dist_res_all <- sum(dist_resid[1:nobs])
  dist_res_newall <- sum(dist_resid_new[1:nobs])
  
  tint_res_all <- sum(tint_resid[1:nobs])
  tint_res_newall <- sum(tint_resid_new[1:nobs])

} # End of the model specification

