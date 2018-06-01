#######################################################################
######################### Nbfragments ######################################
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

Nbfragments.Ms = Nbfragments[voisin =='Ms' & herbi =='E']
Nbfragments.Ju = Nbfragments[voisin =='Ju' & herbi =='E']
data = list(X = Nbfragments.Ms, Y = Nbfragments.Ju, N1 = length(Nbfragments.Ms), N2 = length(Nbfragments.Ju))
boxplot(Nbfragments.Ms, Nbfragments.Ju)

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

boxplot(Nbfragments.Ms, Nbfragments.Ju)
points(1, mean(Ms), col = 'blue')
points(2, mean(Ju), col = 'blue')
points(1, mean(Nbfragments.Ms), col = 'red')
points(2, mean(Nbfragments.Ju), col = 'red')
p.val(Ms)
p.val(Ju)
p.val(Ms - Ju)

mean(Ms)
mean(Ju)
mean(Ms - Ju)

detach(egeria)





########## Jussie ############
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)

Nbfragments.Ms = Nbfragments[voisin =='Ms' & herbi =='E']
Nbfragments.Eg = Nbfragments[voisin =='Eg' & herbi =='E']
data = list(X = Nbfragments.Ms, Y = Nbfragments.Eg, N1 = length(Nbfragments.Ms), N2 = length(Nbfragments.Eg))
boxplot(Nbfragments.Ms, Nbfragments.Eg)

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

boxplot(Nbfragments.Ms, Nbfragments.Eg)

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

Nbfragments.Ju = Nbfragments[voisin =='Ju' & herbi =='E']
Nbfragments.Eg = Nbfragments[voisin =='Eg' & herbi =='E']
data = list(X = Nbfragments.Ju, Y = Nbfragments.Eg, N1 = length(Nbfragments.Ju), N2 = length(Nbfragments.Eg))
boxplot(Nbfragments.Ju, Nbfragments.Eg)

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

boxplot(Nbfragments.Ju, Nbfragments.Eg)
points(rep(1, length(Nbfragments.Ju)), Nbfragments.Ju)
points(rep(2, length(Nbfragments.Eg)), Nbfragments.Eg)
points(1, mean(Nbfragments.Ju, na.rm = TRUE), col = 'red')
points(2, mean(Nbfragments.Eg), col = 'red')
points(1, mean(Ju), col = 'blue')
points(2, mean(Eg), col = 'blue')



p.val(Ju)
p.val(Eg)
p.val(Ju - Eg)

mean(Ju)
mean(Eg)
mean(Ju - Eg)

detach(MSpi)

