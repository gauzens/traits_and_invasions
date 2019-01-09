#######################################################################
######################### RGR ######################################
#######################################################################


setwd('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions')
rm(list = ls())

library(R2jags)
library(car)
source('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/pval.functions.R')


plot.residuals = function(values, design, result){
  predicted = design %*% result
  resids = values - predicted
  resids = (resids - mean(resids, na.rm = TRUE)) / sd(resids, na.rm = TRUE)
  plot(resids ~ predicted, ylab = 'Standardized residuals', xlab = 'Predicted values')
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

model = lm(RGR ~ egeria$herbi * egeria$voisin)
cooks = cooks.distance(model)
hist(cooks, nclass = 20)
# influncial points: greater than 6 times mean of cooks:
cooks[cooks > 6*mean(cooks)]
# corresponding plots: 
plot(cooks)
abline(h=6*mean(cooks))

# test for outliers:
outlierTest(model)

# plot(model)
AIC(model)
AIC(lm(RGR ~ weigths * voisin))
weigths = poidsEcrevisse
# variance of residuals decreases with predicted values. 
# might have forgotten one covariate. bodymass of crayfish might be a good candidate
# as big drayfish might eat more than smaller ones. 
summary(lm(RGR[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass on RGR (when I only consider tratments with crayfish)
# thus, I fit to a simple glm case. No need to add BM as covariate in this case.
# a model RGR ~ biomass * voisin is enough, if I correct heteroscedasticity
# AIC comparison hint in this direction too.
factors = as.numeric(interaction(voisin,herbi))
data = list(y = RGR, X = design, N = dim(egeria)[1], nb.factors = ncol(design), id.factor = factors)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "fixed_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             # n.iter = 500,
                             # n.burnin = 100,
                             DIC = FALSE)
save(model.fit.1, file = '/homes/bg33novu/projects/Lise_ecrevisses/results/Rsave/Egeria_RGR.Rdata')


# intercept = beta[1] = E.JU, mean(RGR[combttt == 'EGJUE'])
# beta[2] = T.JU - E.JU, mean(RGR[combttt == 'EGJUT']) - mean(RGR[combttt == 'EGJUE'])
# beta[3] = E.Ms - E.JU, mean(RGR[combttt == 'EGMSE']) - mean(RGR[combttt == 'EGJUE'])
# beta[4] = T.Ms - beta1 - beta2 - beta3
#         = mean(RGR[combttt == 'EGMST']) - mean(RGR[combttt == 'EGMSE']) - mean(RGR[combttt == 'EGJUT']) + mean(RGR[combttt == 'EGJUE'])
plot(as.mcmc(model.fit.1))
gelman.diag(as.mcmc(model.fit.1))
gelman.plot(as.mcmc(model.fit.1))
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
tapply(RGR, interaction(herbi,voisin), mean, na.rm = TRUE)

predicted = colMeans(means)
boxplot(RGR ~ herbi*voisin)
points(tapply(RGR, interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

# effet interaction:

mu.herbi = 0.5*(T.JU + T.MS) - 0.5*(E.JU + E.MS)
p.herbi = length(mu.herbi[mu.herbi<0]) / length(mu.herbi)
p.val(mu.herbi)

mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.MS + E.MS)



Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)

mu.herbi.egeria = 0.5*(E.JU + E.MS) - 0.5*(T.JU + T.MS)
save(mu.herbi.egeria, file = 'egeria')

detach(egeria)


############# RGR Jussie ############## 
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)

Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(RGR ~ herbi * voisin)

# check model assumptions
model = lm(RGR ~ herbi * voisin)
plot(model)
anova(model)

weigths = poidsEcrevisse
summary(lm(RGR[weigths != 0] ~ weigths[weigths != 0]))
plot(RGR[weigths != 0] ~ weigths[weigths != 0])
# no effect of BM on RGR, go to glm then
AIC(model)
AIC(lm(RGR ~ weigths * voisin))

factors = as.numeric(interaction(voisin,herbi))
data = list(y = RGR, X = design, N = dim(Jussie)[1], nb.factors = ncol(design), id.factor = factors)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "fixed_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)


plot(as.mcmc(model.fit.1))
gelman.diag(as.mcmc(model.fit.1))
gelman.plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))

head(design)


E.EG = results[,1]
T.EG = results[,1] + results[,2]
E.MS = results[,1] + results[,3]
T.MS = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.EG, T.EG, E.MS, T.MS)
names(means) = c('E.EG', 'T.EG', 'E.MS', 'T.MS')

colMeans(means)
tapply(RGR, interaction(herbi,voisin), mean, na.rm = TRUE)


predicted = colMeans(means)
boxplot(RGR ~ herbi*voisin)
points(interaction(herbi,voisin), RGR)
points(tapply(RGR, interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
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

mu.herbi.jussie = 0.5*(E.EG + E.MS) - 0.5*(T.EG + T.MS)
save(mu.herbi.jussie, file = 'jussie')

detach(Jussie)

############## MSpi ############## 
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)

MSpi = tab[tab$sps == 'MSpi', ]
attach(MSpi)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)
boxplot(RGR ~ herbi * voisin)

# one outlier for E.Ju

# check model assumptions, with outlier
model = lm(RGR ~ herbi * voisin)
plot(model)

cooks = cooks.distance(model)
hist(cooks, nclass = 20)
# influncial points: greater than 6 times mean of cooks:
cooks[cooks > 6*mean(cooks)]
# corresponding plots: 
plot(cooks)
abline(h=6*mean(cooks))

# test for outliers:
outlierTest(model)
# plot(model)
# plot(std.res(model) ~ model$fitted.values)
anova(model)
summary(lm(RGR[weigths != 0] ~ weigths[weigths != 0]))
plot(RGR[weigths != 0] ~ weigths[weigths != 0])
# no effect of BM on RGR, go to glm then
AIC(model)
AIC(lm(RGR ~ weigths * voisin))

# check model assumptions, without outlier
model = lm(RGR[-7] ~ herbi[-7] * voisin[-7])
# plot(model)
# plot(std.res(model) ~ model$fitted.values)
# Here residual distribuion is better, but not that good

# then I should do analyses with and without outlier.

### first analyses with the outlier



design = model.matrix( ~  herbi * voisin)
head(design)
factors = as.numeric(interaction(voisin,herbi))
data = list(y = RGR, X = design, N = dim(MSpi)[1], nb.factors = ncol(design), id.factor = factors)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "fixed_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)


plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))

head(design)

E.EG = results[,1]
T.EG = results[,1] + results[,2]
E.JU = results[,1] + results[,3]
T.JU = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.EG, T.EG, E.JU, T.JU)
names(means) = c('E.EG', 'T.EG', 'E.JU', 'T.JU')

colMeans(means)
tapply(RGR, interaction(herbi,voisin), mean, na.rm = TRUE)

predicted = colMeans(means)
boxplot(RGR ~ herbi*voisin)
points(tapply(RGR, interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')



# effet interaction:
mu.herbi = 0.5*(T.JU + T.EG) - 0.5*(E.JU + E.EG)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.EG + E.EG)


Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)



### second, without outlier

herbi2 = herbi[-7]
voisin2 = voisin[-7]
design2 = model.matrix( ~  herbi2 * voisin2)
RGR2 = RGR[-7]
boxplot(RGR2 ~ herbi2 * voisin2)
factors2 = as.numeric(interaction(voisin2,herbi2))
data2 = list(y = RGR2, X = design2, N = length(RGR2), nb.factors = ncol(design), id.factor = factors2)
params = c("beta")
model.fit.2 <- jags.parallel(data = data2,
                             model.file = "fixed_glm.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)


plot(as.mcmc(model.fit.2))
summary2 = summary(as.mcmc(model.fit.2))
summary2$statistics
results2 = as.matrix(as.mcmc(model.fit.2))

head(design)


E.EG = results2[,1]
T.EG = results2[,1] + results2[,2]
E.JU = results2[,1] + results2[,3]
T.JU = results2[,1] + results2[,2]  + results2[,3]  + results2[,4]

means2 = cbind(E.EG, T.EG, E.JU, T.JU)
names(means2) = c('E.EG', 'T.EG', 'E.JU', 'T.JU')

colMeans(means2)
tapply(RGR2, interaction(herbi2,voisin2), mean, na.rm = TRUE)

predicted2 = colMeans(means2)
boxplot(RGR2 ~ herbi2*voisin2)
points(tapply(RGR2, interaction(herbi2,voisin2), mean, na.rm = TRUE), col = 'red')
points(c(predicted2[1], predicted2[2], predicted2[3], predicted2[4]), col = 'blue')



# effet interaction:
mu.herbi2 = 0.5*(T.JU + T.EG) - 0.5*(E.JU + E.EG)
mu.voisin2 = 0.5*(E.JU + T.JU) - 0.5*(T.EG + E.EG)



Interaction2 = results2[,4]
simple.effects(mu.herbi2, mu.voisin2, Interaction2)
pairwise.comps(means2)

mean(mu.herbi2)
mean(mu.voisin2)
mean(Interaction2)

mu.herbi.Mspi = 0.5*(E.JU + E.EG) - 0.5*(T.JU + T.EG)

detach(MSpi)

load('egeria')
load('jussie')

### egeria jussie
mean(mu.herbi.egeria) - mean(mu.herbi.jussie)
p.val(mu.herbi.egeria - mu.herbi.jussie)

### egeria mspi
mean(mu.herbi.egeria) - mean(mu.herbi.Mspi)
p.val(mu.herbi.egeria - mu.herbi.Mspi)


#### jussie mspi
mean(mu.herbi.jussie) - mean(mu.herbi.Mspi)
p.val(mu.herbi.jussie - mu.herbi.Mspi)








