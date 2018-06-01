#######################################################################
######################### SLA ######################################
#######################################################################
setwd('/homes/bg33novu/projects/Lise_ecrevisses/bayesians')
rm(list = ls())
library(R2jags)

plot.residuals = function(values, design, result){
  predicted = design %*% result
  resids = values - predicted
  resids = (resids - mean(resids, na.rm = TRUE)) / sd(resids, na.rm = TRUE)
  plot(resids ~ predicted, ylab = 'Standardized residuals', xlab = 'Predicted values')
}


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
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)

########## EGERIA ############

egeria = tab[tab$sps == 'Egeria', ]
attach(egeria)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(SLA ~ herbi * voisin)
boxplot(log10(SLA) ~ herbi * voisin)
model = lm(log10(SLA) ~ herbi * voisin)
# plot(model)
# plot(std.res(model) ~ model$fitted.values)
# weigths = poidsEcrevisse

model2 = lm(SLA ~ herbi * voisin)
# plot(model2)
# plot(std.res(model2) ~ model2$fitted.values)

# Residuals looks good, simple repeated anova base on log values

summary(lm(LDMC[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(SLA), X = design, N = length(SLA), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             # model.file = 'repeated_anova.txt'
                             model.file = "repeated_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

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
tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE)

gelman.diag(as.mcmc(model.fit.1))
gelman.plot(as.mcmc(model.fit.1))

predicted = colMeans(means)
boxplot(log10(SLA) ~ herbi*voisin)
points(tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

# effet interaction:
mu.herbi = 0.5*(T.JU + T.MS) - 0.5*(E.JU + E.MS)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.MS + E.MS)


Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)

plot.residuals(log10(SLA), design, colMeans(results))

detach(egeria)


##########################3 removing tank nb 56



rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)

########## EGERIA ############

egeria = tab[tab$sps == 'Egeria', ]
outliers = which(egeria$ue == 56)

egeria = egeria[-outliers, ]
attach(egeria)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(SLA ~ herbi * voisin)
boxplot(log10(SLA) ~ herbi * voisin)
model = lm(log10(SLA) ~ herbi * voisin)
# plot(model)
# plot(std.res(model) ~ model$fitted.values)
# weigths = poidsEcrevisse

model2 = lm(SLA ~ herbi * voisin)
# plot(model2)
# plot(std.res(model2) ~ model2$fitted.values)

# Residuals looks good, simple repeated anova base on log values

summary(lm(LDMC[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(SLA), X = design, N = length(SLA), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             # model.file = 'repeated_anova.txt'
                             model.file = "repeated_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

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
tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE)

gelman.diag(as.mcmc(model.fit.1))
gelman.plot(as.mcmc(model.fit.1))

predicted = colMeans(means)
boxplot(log10(SLA) ~ herbi*voisin)
points(tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

# effet interaction:
mu.herbi = 0.5*(T.JU + T.MS) - 0.5*(E.JU + E.MS)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.MS + E.MS)


Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

detach(egeria)



######################################################################3
######################################################################3
######################################################################3




########## Jussie ############

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)

Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(SLA ~ herbi * voisin)
boxplot(log10(SLA) ~ herbi * voisin)
model = lm(log10(SLA) ~ herbi * voisin)
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
# plot(std.res(model) ~ model$fitted.values)
# weigths = poidsEcrevisse


# Residuals looks good, simple repeated anova
summary(lm(SLA[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(SLA), X = design, N = length(SLA), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "repeated_glm.txt",
                             # model.file = "repeated_anova.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

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
tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE)


predicted = colMeans(means)
boxplot(log10(SLA) ~ herbi*voisin)
points(tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')


# effet interaction:
mu.herbi = 0.5*(T.EG + T.MS) - 0.5*(E.EG + E.MS)
mu.voisin = 0.5*(E.EG + T.EG) - 0.5*(T.MS + E.MS)

plot.residuals(log10(SLA), design, colMeans(results))


Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)
detach(Jussie)

# removing otliers, lines 76, 77, 78:

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)
outliers = c(76,77, 78)
Jussie = tab[tab$sps == 'Jussie', ]
Jussie = Jussie[-outliers, ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(SLA ~ herbi * voisin)
boxplot(log10(SLA) ~ herbi * voisin)
model = lm(log10(SLA) ~ herbi * voisin)
# plot(model)
# plot(std.res(model) ~ model$fitted.values)
# weigths = poidsEcrevisse


# Residuals looks good, simple repeated anova
summary(lm(SLA[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(SLA), X = design, N = length(SLA), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "repeated_glm.txt",
                             # model.file = "repeated_anova.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

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
tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE)

predicted = colMeans(means)
boxplot(log10(SLA) ~ herbi*voisin)
points(tapply(log10(SLA), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

# effet interaction:
mu.herbi = 0.5*(T.EG + T.MS) - 0.5*(E.EG + E.MS)
mu.voisin = 0.5*(E.EG + T.EG) - 0.5*(T.MS + E.MS)

Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)

detach(Jussie)
