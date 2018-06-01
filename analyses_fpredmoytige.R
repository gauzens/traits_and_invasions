#######################################################################
######################### Fpredmoytige ######################################
#######################################################################


setwd('/homes/bg33novu/projects/Lise_ecrevisses/bayesians')
rm(list = ls())
library(R2jags)
library(car)


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

Fpredmoytige.Ms = Fpredmoytige[voisin =='Ms' & herbi =='E']/100
Fpredmoytige.Ju = Fpredmoytige[voisin =='Ju' & herbi =='E']/100
# data = list(X = logit(Fpredmoytige.Ms), Y = logit(Fpredmoytige.Ju), N1 = length(Fpredmoytige.Ms), N2 = length(Fpredmoytige.Ju))
data = list(X = Fpredmoytige.Ms, Y = Fpredmoytige.Ju, N1 = length(Fpredmoytige.Ms), N2 = length(Fpredmoytige.Ju))

boxplot(Fpredmoytige.Ms, Fpredmoytige.Ju)
points(rep(1, length(Fpredmoytige.Ms)),Fpredmoytige.Ms)
points(rep(2, length(Fpredmoytige.Ju)),Fpredmoytige.Ju)
boxplot(logit(Fpredmoytige.Ms), logit(Fpredmoytige.Ju))
points(rep(1, length(Fpredmoytige.Ms)),logit(Fpredmoytige.Ms))
points(rep(2, length(Fpredmoytige.Ju)),logit(Fpredmoytige.Ju))

params = c('alpha', 'beta', 'sigma.alpha', 'sigma.beta')
model.fit.1 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 200000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))

Ms = results[,1]
Ju = results[,2]

gelman.diag(as.mcmc(model.fit.1))
gelman.plot(as.mcmc(model.fit.1))
# 
# boxplot(Fpredmoytige.Ms, Fpredmoytige.Ju)
# boxplot(logit(Fpredmoytige.Ms), logit(Fpredmoytige.Ju))
# points(rep(1, length(Fpredmoytige.Ms)),logit(Fpredmoytige.Ms))
# points(rep(2, length(Fpredmoytige.Ju)),logit(Fpredmoytige.Ju))
# 
# points(c(mean(Fpredmoytige.Ms, na.rm = TRUE), mean(Fpredmoytige.Ju, na.rm = TRUE)), col = 'red')
# points(c(mean(Ms), mean(Ju)), col = 'blue')
# 
# rbind(logit(Fpredmoytige.Ms), Fpredmoytige.Ms)
# rbind(logit(Fpredmoytige.Ju), Fpredmoytige.Ju)

hist(Ms - Ju)
hist(Ms)
hist(Ju)

p.val(Ms)
p.val(Ju)
p.val(Ms - Ju)


detach(egeria)

###########333 try to do an anova

egeria = tab[tab$sps == 'Egeria', ]
attach(egeria)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

factors = as.numeric(interaction(voisin,herbi))

model = lm(Fpredmoytige ~ egeria$herbi * egeria$voisin)

data = list(y = Fpredmoytige, X = design, N = dim(egeria)[1], nb.factors = ncol(design), id.factor = factors)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "fixed_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)



# intercept = beta[1] = E.JU, mean(RGR[combttt == 'EGJUE'])
# beta[2] = T.JU - E.JU, mean(RGR[combttt == 'EGJUT']) - mean(RGR[combttt == 'EGJUE'])
# beta[3] = E.Ms - E.JU, mean(RGR[combttt == 'EGMSE']) - mean(RGR[combttt == 'EGJUE'])
# beta[4] = T.Ms - beta1 - beta2 - beta3
#         = mean(RGR[combttt == 'EGMST']) - mean(RGR[combttt == 'EGMSE']) - mean(RGR[combttt == 'EGJUT']) + mean(RGR[combttt == 'EGJUE'])
plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))
head(results)

E.JU = results[,1]
T.JU = results[,1] + results[,2]
E.MS = results[,1] + results[,3]
T.MS = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.JU, T.JU, E.MS, T.MS)
names(means) = c('E.JU', 'T.JU', 'E.MS', 'T.MS')
colMeans(means)
tapply(Fpredmoytige, interaction(herbi,voisin), mean, na.rm = TRUE)



predicted = colMeans(means)
boxplot(Fpredmoytige ~ herbi*voisin)
points(tapply(Fpredmoytige, interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

mu.herbi = 0.5*(T.JU + T.MS) - 0.5*(E.JU + E.MS)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.MS + E.MS)
Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)
# effet interaction:

detach(egeria)


########## Jussie ############
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)

Fpredmoytige.Ms = Fpredmoytige[voisin =='Ms' & herbi =='E']/100
Fpredmoytige.Eg = Fpredmoytige[voisin =='Eg' & herbi =='E']/100
data = list(X = Fpredmoytige.Ms, Y = Fpredmoytige.Eg, N1 = length(Fpredmoytige.Ms), N2 = length(Fpredmoytige.Eg))
boxplot(Fpredmoytige.Ms, Fpredmoytige.Eg)

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

##################################### anova

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

factors = as.numeric(interaction(voisin,herbi))

model = lm(Fpredmoytige ~ Jussie$herbi * Jussie$voisin)

data = list(y = Fpredmoytige, X = design, N = length(Fpredmoytige), nb.factors = ncol(design), id.factor = factors)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "fixed_glm_0var.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 50000,
                             n.burnin = 10000,
                             DIC = FALSE)



# intercept = beta[1] = E.JU, mean(RGR[combttt == 'EGJUE'])
# beta[2] = T.JU - E.JU, mean(RGR[combttt == 'EGJUT']) - mean(RGR[combttt == 'EGJUE'])
# beta[3] = E.Ms - E.JU, mean(RGR[combttt == 'EGMSE']) - mean(RGR[combttt == 'EGJUE'])
# beta[4] = T.Ms - beta1 - beta2 - beta3
#         = mean(RGR[combttt == 'EGMST']) - mean(RGR[combttt == 'EGMSE']) - mean(RGR[combttt == 'EGJUT']) + mean(RGR[combttt == 'EGJUE'])
plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))
head(results)

E.JU = results[,1]
T.JU = results[,1] + results[,2]
E.MS = results[,1] + results[,3]
T.MS = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.JU, T.JU, E.MS, T.MS)
names(means) = c('E.JU', 'T.JU', 'E.MS', 'T.MS')
colMeans(means)
tapply(Fpredmoytige, interaction(herbi,voisin), mean, na.rm = TRUE)

predicted = colMeans(means)
boxplot(Fpredmoytige ~ herbi*voisin)
points(tapply(Fpredmoytige, interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

mu.herbi = 0.5*(T.JU + T.MS) - 0.5*(E.JU + E.MS)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.MS + E.MS)
Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)
# effet interaction:

detach(Jussie)




########## MSpi ############
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
MSpi = tab[tab$sps == 'MSpi', ]
attach(MSpi)

Fpredmoytige.Ju = Fpredmoytige[voisin =='Ju' & herbi =='E']/100
Fpredmoytige.Eg = Fpredmoytige[voisin =='Eg' & herbi =='E']/100
data = list(X = Fpredmoytige.Ju, Y = Fpredmoytige.Eg, N1 = length(Fpredmoytige.Ju), N2 = length(Fpredmoytige.Eg))
boxplot(Fpredmoytige.Ju, Fpredmoytige.Eg)
points(rep(1, length(Fpredmoytige.Ju)), Fpredmoytige.Ju)
points(rep(2, length(Fpredmoytige.Eg)), Fpredmoytige.Eg)

paraJu = c('alpha', 'beta', 'sigma.alpha', 'sigma.beta')
model.fit.1 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = paraJu,
                             n.chains = 4,
                             n.iter = 200000,
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

boxplot(Fpredmoytige.Ju, Fpredmoytige.Eg)
# points(rep(1, length(Fpredmoytige.Ju)), Fpredmoytige.Ju)
# points(rep(2, length(Fpredmoytige.Eg)), Fpredmoytige.Eg)
points(1, mean(Fpredmoytige.Ju, na.rm = TRUE), col = 'red')
points(2, mean(Fpredmoytige.Eg), col = 'red')
points(1, mean(Ju), col = 'blue')
points(2, mean(Eg), col = 'blue')
p.val(Ju)
p.val(Eg)
p.val(Ju - Eg)

mean(Ju)
mean(Eg)
p.val(Ju - Eg)

detach(MSpi)




# anova ace espece et herbi en facteur



rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
tab2 = tab[tab$presence == 'Herbi',]
attach(tab2)

sp = droplevels(sps)
herbi = droplevels(herbi)
design = model.matrix( ~  sp * voisin)
colSums(design)
factors = as.numeric(interaction(sps,voisin))
model = lm(RGR ~ tab$sp * tab$voisin)


data = list()

data = list(y = Fpredmoytige, X = design, N = length(Fpredmoytige), nb.factors = ncol(design), id.factor = factors)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "fixed_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 50000,
                             n.burnin = 10000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))
head(results)

p.val(results[,8])
 
