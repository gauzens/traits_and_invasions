#######################################################################
######################### Freebiomass ######################################
#######################################################################


setwd('/homes/bg33novu/projects/Lise_ecrevisses/bayesians')
rm(list = ls())
library(R2jags)


std.res = function(model){
  return((model$residuals - mean(model$residuals))/sd(model$residuals))
}

p.val = function(vec){
  mean.val = mean(vec)
  if (mean.val > 0){
    p.value = length(vec[vec<0])/length(vec)
  }else{
    p.value = length(vec[vec>0])/length(vec)
  }
  return(p.value)
}


rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)


########## EGERIA ############

egeria = tab[tab$sps == 'Egeria', ]
attach(egeria)

Freebiomass.Ms = Freebiomass[voisin =='Ms' & herbi =='E']
Freebiomass.Ju = Freebiomass[voisin =='Ju' & herbi =='E']
data = list(X = Freebiomass.Ms, Y = Freebiomass.Ju, N1 = length(Freebiomass.Ms), N2 = length(Freebiomass.Ju))
boxplot(Freebiomass.Ms, Freebiomass.Ju)

params = c('alpha', 'beta', 'sigma.alpha', 'sigma.beta')
model.fit.1 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))

Ms = results[,1]
Ju = results[,2]

hist(Ms - Ju)
hist(Ms)
hist(Ju)

p.val(Ms)
p.val(Ju)
p.val(Ms - Ju)

mean(Ms)
mean(Ju)
mean(Ms - Ju)

gelman.diag(as.mcmc(model.fit.1))
gelman.plot(as.mcmc(model.fit.1))
detach(egeria)





########## Jussie ############
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)

Freebiomass.Ms = Freebiomass[voisin =='Ms' & herbi =='E']
Freebiomass.Eg = Freebiomass[voisin =='Eg' & herbi =='E']
data = list(X = Freebiomass.Ms, Y = Freebiomass.Eg, N1 = length(Freebiomass.Ms), N2 = length(Freebiomass.Eg))
boxplot(Freebiomass.Ms, Freebiomass.Eg)

params = c('alpha', 'beta', 'sigma.alpha', 'sigma.beta')
model.fit.1 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))

Ms = results[,1]
Eg = results[,2]

hist(Ms - Eg)
hist(Ms)
hist(Eg)

p.val(Ms)
p.val(Eg)
p.val(Ms - Eg)

mean(Ms)
mean(Eg)
mean(Ms - Eg)

detach(Jussie)





########## MSpi ############
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
MSpi = tab[tab$sps == 'MSpi', ]
attach(MSpi)

Freebiomass.Ju = Freebiomass[voisin =='Ju' & herbi =='E']
Freebiomass.Eg = Freebiomass[voisin =='Eg' & herbi =='E']
data = list(X = Freebiomass.Ju, Y = Freebiomass.Eg, N1 = length(Freebiomass.Ju), N2 = length(Freebiomass.Eg))
boxplot(Freebiomass.Ju, Freebiomass.Eg)
boxplot(log(Freebiomass.Ju+1), log(Freebiomass.Eg+1))


paraJu = c('alpha', 'beta', 'sigma.alpha', 'sigma.beta')
model.fit.1 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = paraJu,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))

Ju = results[,1]
Eg = results[,2]

hist(Ju - Eg)
hist(Ju)
hist(Eg)

p.val(Ju)
p.val(Eg)
p.val(Ju - Eg)

mean(Ju)
mean(Eg)
mean(Ju - Eg)

detach(MSpi)

