
# 1 Anova, to be used for: RGR, DMCBiomass
cat("
    model
    {
      # priors
      alpha ~ dnorm(0.006, 1.0E-6)
      beta ~ dnorm(0.006, 1.0E-6)
      # residuals, group specific
      sigma.alpha ~ dunif(0,100) # (you may want to use a more proper gamma prior)
      tau.alpha <- 1/(sigma.alpha*sigma.alpha)
      sigma.beta ~ dunif(0,100) # (you may want to use a more proper gamma prior)
      tau.beta <- 1/(sigma.beta*sigma.beta)

      # likelihood
      for(i in 1:N1)
      {
        X[i] ~ dnorm(alpha, tau.alpha)
      }
      for(i in 1:N2)
      {
        mean[i] = alpha
        Y[i] ~ dnorm(beta, tau.beta)
      }
    
    }
    
    ", file="/homes/bg33novu/projects/Lise_ecrevisses/bayesians/Welch_t.test.txt")

