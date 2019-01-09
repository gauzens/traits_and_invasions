#######################################################################
######################### LDMC ######################################
#######################################################################
setwd('/homes/bg33novu/projects/Lise_ecrevisses/bayesians')
library(R2jags)
library(car)

rm(list = ls())


source('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/pval.functions.R')
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)

########## EGERIA ############

egeria = tab[tab$sps == 'Egeria', ]
attach(egeria)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(LDMC ~ herbi * voisin)
boxplot(log10(LDMC) ~ herbi * voisin)
library(nlme)
# model = lme(log10(LDMC) ~ herbi * voisin, random = ~1|ue, data = egeria, na.action = na.omit)
model = lm(log10(LDMC) ~ herbi * voisin)

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
summary(lm(LDMC[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(LDMC), X = design, N = dim(egeria)[1], nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             # model.file = 'repeated_anova.txt', 
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
tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE)

# effet interaction:
mu.herbi = 0.5*(T.JU + T.MS) - 0.5*(E.JU + E.MS)

mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.MS + E.MS)


boxplot(log10(LDMC) ~ voisin*herbi)
predicted = colMeans(means)
boxplot(log10(LDMC) ~ herbi*voisin)
points(tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)


gelman.diag(as.mcmc(model.fit.1))
gelman.plot(as.mcmc(model.fit.1))

plot.residuals(log10(LDMC), design, colMeans(results))

############## without outliers
# outliers = c(16, 21, 65)
outliers = 16
voisin2 = droplevels(voisin[-outliers])
herbi2 = droplevels(herbi[-outliers])
design2 = model.matrix( ~  herbi2 * voisin2)
colSums(design2)
LDMC2 = LDMC[-outliers]
model2 = lm(log10(LDMC2) ~ herbi2 * voisin2)
# plot(model2)


factors2 = as.numeric(interaction(voisin2,herbi2))
tanks2 = paste('t', ue[-outliers], sep = '')
tanks2 = as.numeric(as.factor(as.character(ue[-outliers])))
nb_tanks2 = max(tanks2[-outliers])
data2 = list(y = log10(LDMC2), X = design2, N =length(LDMC2), nb.factors = ncol(design2), id.factor = factors2, nb_tanks = nb_tanks2, id = tanks2)
params = c("beta")
model.fit.2 <- jags.parallel(data = data2,
                             # model.file = "repeated_anova.txt",
                             model.file = "repeated_glm.txt", 
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)


plot(as.mcmc(model.fit.2))
summary2 = summary(as.mcmc(model.fit.2))
summary2$statistics
results2 = as.matrix(as.mcmc(model.fit.2))
head(results2)

E.JU2 = results2[,1]
T.JU2 = results2[,1] + results2[,2]
E.MS2 = results2[,1] + results2[,3]
T.MS2 = results2[,1] + results2[,2]  + results2[,3]  + results2[,4]

means2 = cbind(E.JU2, T.JU2, E.MS2, T.MS2)
names(means2) = c('E.JU', 'T.JU', 'E.MS', 'T.MS')
colMeans(means2)
tapply(log10(LDMC2), interaction(herbi2,voisin2), mean, na.rm = TRUE)

# effet interaction:
mu.herbi2 = 0.5*(T.JU2 + T.MS2) - 0.5*(E.JU2 + E.MS2)
mu.voisin2 = 0.5*(E.JU2 + T.JU2) - 0.5*(T.MS2 + E.MS2)


gelman.diag(as.mcmc(model.fit.2))
gelman.plot(as.mcmc(model.fit.2))

predicted2 = colMeans(means2)
boxplot(log10(LDMC2) ~ herbi2*voisin2)
points(log10(LDMC2) ~ interaction(herbi2,voisin2))
points(tapply(log10(LDMC2), interaction(herbi2,voisin2), mean, na.rm = TRUE), col = 'red')
points(c(predicted2[1], predicted2[2], predicted2[3], predicted2[4]), col = 'blue')


Interaction2 = results2[,4]
simple.effects(mu.herbi2, mu.voisin2, Interaction2)
pairwise.comps(means2)

mean(mu.herbi2)
mean(mu.voisin2)
mean(Interaction2)


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

boxplot(LDMC ~ herbi * voisin)
boxplot(log10(LDMC) ~ herbi * voisin)
model = lm(log10(LDMC) ~ herbi * voisin)
# plot(model)
# plot(std.res(model) ~ model$fitted.values)
# weigths = poidsEcrevisse

cooks = cooks.distance(model)
hist(cooks, nclass = 20)
# influncial points: greater than 6 times mean of cooks:
cooks[cooks > 6*mean(cooks)]
# corresponding plots: 
plot(cooks)
abline(h=6*mean(cooks))


# test for outliers:
outlierTest(model)

summary(lm(LDMC[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(LDMC), X = design, N = length(LDMC), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "repeated_glm.txt",
                             # model.file = "repeated_anova.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 200000,
                             n.burnin = 100000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1), start = 20000)
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
tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE)
boxplot(log10(LDMC) ~ voisin*ttt)
points(tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points()

(colMeans(means) - tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE))/colMeans(means)

points(c(mean(E.EG), mean(E.MS), mean(T.EG), mean(T.MS)), col = 'red')


hist(T.EG)
points(mean(log10(LDMC)[herbi == 'T' & voisin == 'Egeria']), 0, col = 'red')
mean(T.EG)
mean(log10(LDMC)[herbi == 'T' & voisin == 'Egeria'])
# effet interaction:
mu.herbi = 0.5*(T.EG + T.MS) - 0.5*(E.EG + E.MS)
mu.voisin = 0.5*(E.EG + T.EG) - 0.5*(T.MS + E.MS)

as.mcmc(model.fit.1)[[3]]
boxplot(log10(LDMC) ~ voisin*ttt)
tapply(log10(LDMC), interaction(voisin,herbi), mean, na.rm = TRUE)

predicted = colMeans(means)
boxplot(log10(LDMC) ~ voisin*ttt)
points(tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[3], predicted[2], predicted[4]), col = 'blue')

Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)
detach(Jussie)

# if I remove ue56 pluis one big oultier line 64:
# their LDMC values was 16.4, 15.6 and 19
# in Egeria, Mspi and Egeria

detach(Jussie)

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)

Jussie = tab[tab$sps == 'Jussie', ]
ue56 = which(Jussie$ue == 56)
Jussie = Jussie[-c(ue56, 64), ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(LDMC ~ herbi * voisin)

boxplot(tab[tab$sps == 'Jussie', ]$LDMC ~ tab[tab$sps == 'Jussie', ]$herbi * tab[tab$sps == 'Jussie', ]$voisin)

boxplot(log10(LDMC) ~ herbi * voisin)
aa = tab[tab$sps == 'Jussie', ]
boxplot(log10(aa$LDMC) ~ droplevels(aa$herbi) * droplevels(aa$voisin))
model = lm(log10(LDMC) ~ herbi * voisin)
model2 = lm(LDMC ~ herbi * voisin)
plot(model)
plot(model2)
plot(std.res(model) ~ model$fitted.values)
plot(std.res(model2) ~ model2$fitted.values)
# weigths = poidsEcrevisse


# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(LDMC), X = design, N = length(LDMC), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "repeated_glm.txt", #repeated_anova.txt
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 200000,
                             n.burnin = 100000,
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
tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE)

# effet interaction:
mu.herbi = 0.5*(T.EG + T.MS) - 0.5*(E.EG + E.MS)
mu.voisin = 0.5*(E.EG + T.EG) - 0.5*(T.MS + E.MS)

boxplot(LDMC ~ voisin*herbi)
Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)


detach(Jussie)








########## Mspi ############

rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)


Mspi = tab[tab$sps == 'Mspi', ]
attach(Mspi)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)
colSums(design)

boxplot(LDMC ~ herbi * voisin)
boxplot(log10(LDMC) ~ herbi * voisin)


boxplot(log10(LDMC) ~ voisin*ttt)
points(tapply(log10(LDMC), interaction(voisin,ttt), mean, na.rm = TRUE), col = 'red')
points(interaction(voisin,ttt), log10(LDMC))
model = lm(log10(LDMC) ~ voisin*ttt)
# plot(model)
# plot(std.res(model) ~ model$fitted.values)
# weigths = poidsEcrevisse
# model = lm(LDMC ~ herbi * voisin)
# plot(model)

cooks = cooks.distance(model)
hist(cooks, nclass = 20)
# influncial points: greater than 6 times mean of cooks:
cooks[cooks > 6*mean(cooks)]
# corresponding plots: 
plot(cooks)
abline(h=6*mean(cooks))

# test for outliers:
library(car)
outlierTest(model)


# Residuals looks good, simple repeated anova
summary(lm(LDMC[weigths != 0] ~ weigths[weigths != 0]))
# no effect of crayfisch bodymass 

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(LDMC), X = design, N = length(LDMC), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.fit.1 <- jags.parallel(data = data,
                             model.file = "repeated_glm.txt", #repeated_anova.txt
                             # model.file = "repeated_anova.txt",
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 200000,
                             n.burnin = 100000,
                             DIC = FALSE)


plot(as.mcmc(model.fit.1))
summary = summary(as.mcmc(model.fit.1))
summary$statistics
results = as.matrix(as.mcmc(model.fit.1))
head(results)

plot.residuals(log10(LDMC), design, colMeans(results))
  


E.EG = results[,1]
T.EG = results[,1] + results[,2]
E.JU = results[,1] + results[,3]
T.JU = results[,1] + results[,2]  + results[,3]  + results[,4]

means = cbind(E.EG, T.EG, E.JU, T.JU)
names(means) = c('E.EG', 'T.EG', 'E.JU', 'T.JU')
predicted = colMeans(means)
colMeans(means)
tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE)

boxplot(log10(LDMC) ~ ttt*voisin)
points(tapply(log10(LDMC), interaction(herbi,voisin), mean, na.rm = TRUE), col = 'red')
points(c(predicted[1], predicted[2], predicted[3], predicted[4]), col = 'blue')

# effet interaction:
mu.herbi = 0.5*(T.JU + T.EG) - 0.5*(E.JU + E.EG)
mu.voisin = 0.5*(E.JU + T.JU) - 0.5*(T.EG + E.EG)

Interaction = results[,4]
simple.effects(mu.herbi, mu.voisin, Interaction)
pairwise.comps(means)

mean(mu.herbi)
mean(mu.voisin)
mean(Interaction)


detach(Mspi)

