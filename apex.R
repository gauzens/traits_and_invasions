#######################################################################
######################### RGR ######################################
#######################################################################


setwd('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions')
rm(list = ls())

library(R2jags)
library(car)
library(MASS)
source('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/pval.functions.R')

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

plot.residuals = function(values, design, result){
  predicted = design %*% result
  resids = values - predicted
  resids = (resids - mean(resids, na.rm = TRUE)) / sd(resids, na.rm = TRUE)
  plot(resids ~ predicted, ylab = 'Standardized residuals', xlab = 'Predicted values')
}

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/apex.csv', header = T)
tab$apex = tab$apex / 100
tab$voisin = tab$spscomp
tab$herbi = tab$ttt
tab$apex.logit = logit(tab$apex)

########## EGERIA ############

egeria = tab[tab$sps == 'Egeria', ]
attach(egeria)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)
model = lm(apex.logit ~ egeria$herbi * egeria$voisin)
plot(model)

apex.Ms = apex[voisin =='Ms' & herbi =='E']
apex.Ju = apex[voisin =='Ju' & herbi =='E']

data = list(X = apex.Ms, Y = apex.Ju, N1 = length(apex.Ms), N2 = length(apex.Ju))
boxplot(apex.Ms, apex.Ju)

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


boxplot(apex.Ms, apex.Ju)
points(1, mean(Ms), col = 'blue')
points(2, mean(Ju), col = 'blue')
points(1, mean(apex.Ms), col = 'red')
points(2, mean(apex.Ju), col = 'red')
p.val(Ms-1)
p.val(Ju-1)
p.val(Ms - Ju)

mean(Ms)
mean(Ju)
mean(Ms - Ju)


detach(egeria)





########## Jussie ############
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/apex.csv', header = T)
tab$apex = tab$apex / 100
tab$voisin = tab$spscomp
tab$herbi = tab$ttt
tab$apex.logit = logit(tab$apex)


Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)

apex.Ms = apex[voisin =='Ms' & herbi =='E']
apex.Eg = apex[voisin =='Eg' & herbi =='E']
data = list(X = apex.Ms, Y = apex.Eg, N1 = length(apex.Ms), N2 = length(apex.Eg))
boxplot(apex.Ms, apex.Eg)

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

boxplot(apex.Ms, apex.Eg)

p.val(Ms)
p.val(Eg)
p.val(Ms - Eg)

mean(Ms)
mean(Eg)
mean(Ms - Eg)

detach(Jussie)



########## MSpi ############
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/apex.csv', header = T)
tab$apex = tab$apex / 100
tab$voisin = tab$spscomp
tab$herbi = tab$ttt
tab$apex.logit = logit(tab$apex)


MSpi = tab[tab$sps == 'MSpi', ]
attach(MSpi)

apex.Ju = apex[voisin =='Ju' & herbi =='E']
apex.Eg = apex[voisin =='Eg' & herbi =='E']
data = list(X = apex.Ju, Y = apex.Eg, N1 = length(apex.Ju), N2 = length(apex.Eg))
boxplot(apex.Ju, apex.Eg)

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

Ju = results[,1]
Eg = results[,2]

hist(Ju - Eg)
hist(Ju)
hist(Eg)

boxplot(apex.Ju, apex.Eg)

p.val(Ju-1)
p.val(Eg-1)
p.val(Ju - Eg)

mean(Ju)
mean(Eg)
mean(Ms - Eg)

detach(MSpi)


