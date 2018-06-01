#######################################################################
######################### DMCbiom ######################################
#######################################################################



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

simple.effects = function(a, b, c){
  cat('---------- Simple effects: ---------\n')
  cat(deparse(substitute(a)), ': ', p.val(a), '\n')
  cat(deparse(substitute(b)), ': ', p.val(b), '\n')
  cat(deparse(substitute(c)), ': ', p.val(c), '\n')
  
  cat('------------------------------------\n')
}

pairwise.comps = function(means){
  cat('---------- Pairwise comparisons: ---------\n')
  for (i in 1:dim(means)[2]){
    if (i<dim(means)[2]){
      for (j in (i+1):dim(means)[2]){
        cat(names(means)[i], ' - ', names(means)[j],': ', format(p.val(means[,i] - means[,j]), digits = 6), '\n')
      }
    }
  }
  cat('------------------------------------------\n')
}






rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)

########## EGERIA ############

egeria = tab[tab$sps == 'Egeria', ]
attach(egeria)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(DMCbiom ~ herbi * voisin)
model = lm(DMCbiom ~ herbi * voisin)
# plot(model)
weigths = poidsEcrevisse

# variance: T.Ju looks lower than other groups
summary(lm(DMCbiom[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
data = list(y = DMCbiom, X = design, N = dim(egeria)[1], nb.factors = ncol(design), id.factor = factors)
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
tapply(DMCbiom, interaction(herbi,voisin), mean, na.rm = TRUE)

# effet interaction:
mu.herbi = 0.5*(T.JU + T.MS) - 0.5*(E.JU + E.MS)
hist(mu.herbi)
p.herbi = length(mu.herbi[mu.herbi<0]) / length(mu.herbi)
p.val(mu.herbi)

mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.MS + E.MS)
hist(mu.voisin)

boxplot(DMCbiom ~ voisin*herbi)
predicted = colMeans(means)
boxplot(DMCbiom ~ herbi*voisin)
points(tapply(DMCbiom, interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)


detach(egeria)



########## Jussie ############


rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)
boxplot(DMCbiom ~ herbi * voisin)

model = lm(DMCbiom ~ herbi * voisin)
# plot(model)
weigths = poidsEcrevisse

# variance looks not so great, maybe some effects of herbivory on resiual variance.
# would use glm
summary(lm(DMCbiom[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
data = list(y = DMCbiom, X = design, N = length(DMCbiom), nb.factors = ncol(design), id.factor = factors)
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

E.EG = results[,1]
T.EG = results[,1] + results[,2]
E.MS = results[,1] + results[,3]
T.MS = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.EG, T.EG, E.MS, T.MS)
names(means) = c('E.EG', 'T.EG', 'E.MS', 'T.MS')

colMeans(means)
tapply(LDMC, interaction(herbi,voisin), mean, na.rm = TRUE)

# effet interaction:
mu.herbi = 0.5*(T.EG + T.MS) - 0.5*(E.EG + E.MS)
mu.voisin = 0.5*(E.EG + T.EG) - 0.5*(T.MS + E.MS)


boxplot(DMCbiom ~ voisin*herbi)
Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)


detach(Jussie)


################ MSpi ##############################

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
MSpi = tab[tab$sps == 'MSpi', ]
attach(MSpi)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)
boxplot(log10(DMCbiom) ~ herbi * voisin)
boxplot(DMCbiom ~ herbi * voisin)

model = lm(DMCbiom ~ herbi * voisin)
# plot(model)

cooks = cooks.distance(model)
hist(cooks, nclass = 20)
# influncial points: greater than 6 times mean of cooks:
cooks[cooks > 6*mean(cooks)]
# corresponding plots: 
plot(cooks)
abline(h=6*mean(cooks))

# test for outliers:
outlierTest(model)
weigths = poidsEcrevisse

# variance looks not so great, maybe some effects of neighboorhood on resiual variance.
# would use glm
summary(lm(DMCbiom[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
data = list(y = DMCbiom, X = design, N = length(DMCbiom), nb.factors = ncol(design), id.factor = factors)
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

E.EG = results[,1]
T.EG = results[,1] + results[,2]
E.JU = results[,1] + results[,3]
T.JU = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.EG, T.EG, E.JU, T.JU)
names(means) = c('E.EG', 'T.EG', 'E.JU', 'T.JU')

colMeans(means)
tapply(DMCbiom, interaction(herbi,voisin), mean, na.rm = TRUE)

# effet interaction:
mu.herbi = 0.5*(T.JU + T.EG) - 0.5*(E.JU + E.EG)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.EG + E.EG)

Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)


detach(MSpi)

################ removing outliers

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
outliers = c(4, 10, 11)
MSpi = tab[tab$sps == 'MSpi', ]
MSpi = MSpi[-outliers, ]

attach(MSpi)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)
boxplot(DMCbiom ~ herbi * voisin)

model = lm(DMCbiom ~ herbi * voisin)
# plot(model)
weigths = poidsEcrevisse

# variance looks not so great, maybe some effects of neighboorhood on resiual variance.
# would use glm
summary(lm(DMCbiom[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
data = list(y = DMCbiom, X = design, N = length(DMCbiom), nb.factors = ncol(design), id.factor = factors)
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

E.EG = results[,1]
T.EG = results[,1] + results[,2]
E.JU = results[,1] + results[,3]
T.JU = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.EG, T.EG, E.JU, T.JU)
names(means) = c('E.EG', 'T.EG', 'E.JU', 'T.JU')

colMeans(means)
tapply(DMCbiom, interaction(herbi,voisin), mean, na.rm = TRUE)

# effet interaction:
mu.herbi = 0.5*(T.JU + T.EG) - 0.5*(E.JU + E.EG)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.EG + E.EG)

Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)

detach(MSpi)


