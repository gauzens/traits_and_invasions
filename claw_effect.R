#### comparaisons pinces sans pinces ########
rm(list = ls())
library(R2jags)
setwd('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/')
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



tab.esp = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/Traits_ESP.csv', header = T, sep = ',')
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv')

sink('claws_comparisons.txt')
cat('################################################\n')
cat('######### effect of claw removal ###############\n')
cat('################################################\n\n\n\n\n')
sink()

##################################################
##################### Egeria #####################
##################################################

sink('claws_comparisons.txt', append = TRUE)
cat('####################\n')
cat('###### Egeria ######\n')
cat('####################\n\n')
sink()
##########333 nb fragments ############
rm(list = setdiff(ls(), lsf.str()))
tab.esp = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/Traits_ESP.csv', header = T, sep = ',')
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv')
with.claw.frags = tab$Nbfragments[tab$herbi == 'E' & tab$sps == 'Egeria']
without.claw.frags = tab.esp$Nbfragments[tab.esp$sps == 'Egeria']
base = tab$Nbfragments[tab$herbi == 'T' & tab$sps == 'Egeria']

## pinces vs sans pinces ###
data = list(X = with.claw.frags, Y = without.claw.frags, N1 = length(with.claw.frags), N2 = length(without.claw.frags))


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

distrib.with.claw = results[,1]
distrib.without.claw = results[,2]

## pvalue with claws, without claws
p.val(results[,1] - results[,2])


data = list(X = without.claw.frags, Y = base, N1 = length(with.claw.frags), N2 = length(base))
model.fit.2 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.2))
summary = summary(as.mcmc(model.fit.2))
summary$statistics
results2 = as.matrix(as.mcmc(model.fit.2))

## pvalue without claws, without crayfish
p.val(results2[,1] - results2[,2])


data = list(X = with.claw.frags, Y = base, N1 = length(with.claw.frags), N2 = length(without.claw.frags))
model.fit.3 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.3))
summary = summary(as.mcmc(model.fit.3))
summary$statistics
results3 = as.matrix(as.mcmc(model.fit.3))

sink('claws_comparisons.txt', append = TRUE)
cat('###nb.fragments\n\n')

cat('pvalue with claws vs without claws\n')
p.val(results[,1] - results[,2])
cat('\n')

cat('pvalue without claws vs temoin\n')
p.val(results2[,1] - results2[,2])
cat('\n')

cat('pvalue with claws vs temoin\n')
p.val(results3[,1] - results3[,2])
cat('\n')

sink()


##### LDMC ###############
rm(list = setdiff(ls(), lsf.str()))
tab.esp = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/Traits_ESP_LDMC.csv', header = T)
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)
with.claw.LDMC = tab$LDMC[tab$herbi == 'E' & tab$sps == 'Egeria']
without.claw.LDMC = tab.esp$LDMC[tab.esp$sps == 'Egeria']
base = tab$Nbfragments[tab$herbi == 'T' & tab$sps == 'Egeria']

## pinces vs sans pinces ###
data = list(X = with.claw.LDMC, Y = without.claw.LDMC, N1 = length(with.claw.LDMC), N2 = length(without.claw.LDMC))


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

distrib.with.claw = results[,1]
distrib.without.claw = results[,2]

## pvalue with claws, without claws
p.val(results[,1] - results[,2])


data = list(X = without.claw.LDMC, Y = base, N1 = length(without.claw.LDMC), N2 = length(base))
model.fit.2 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.2))
summary = summary(as.mcmc(model.fit.2))
summary$statistics
results2 = as.matrix(as.mcmc(model.fit.2))

## pvalue without claws, without crayfish 
p.val(results2[,1] - results2[,2])


data = list(X = with.claw.LDMC, Y = base, N1 = length(with.claw.LDMC), N2 = length(base))
model.fit.3 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.3))
summary = summary(as.mcmc(model.fit.3))
summary$statistics
results3 = as.matrix(as.mcmc(model.fit.3))


sink('claws_comparisons.txt', append = TRUE)
cat('###LDMC\n\n')

cat('pvalue with claws vs without claws\n')
p.val(results[,1] - results[,2])
cat('\n')

cat('pvalue without claws vs temoin\n')
p.val(results2[,1] - results2[,2])
cat('\n')

cat('pvalue with claws vs temoin\n')
p.val(results3[,1] - results3[,2])
cat('\n')

sink()


########################## jussie #####################################



sink('claws_comparisons.txt', append = TRUE)
cat('####################\n')
cat('###### Jussie ######\n')
cat('####################\n\n')
sink()
##########333 nb fragments ############
rm(list = setdiff(ls(), lsf.str()))
tab.esp = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/Traits_ESP.csv', header = T, sep = ',')
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv')
with.claw.frags = tab$Nbfragments[tab$herbi == 'E' & tab$sps == 'Jussie']
without.claw.frags = tab.esp$Nbfragments[tab.esp$sps == 'Jussie']
base = tab$Nbfragments[tab$herbi == 'T' & tab$sps == 'Jussie']

## pinces vs sans pinces ###
data = list(X = with.claw.frags, Y = without.claw.frags, N1 = length(with.claw.frags), N2 = length(without.claw.frags))


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

distrib.with.claw = results[,1]
distrib.without.claw = results[,2]

## pvalue with claws, without claws
p.val(results[,1] - results[,2])


data = list(X = without.claw.frags, Y = base, N1 = length(with.claw.frags), N2 = length(without.claw.frags))
model.fit.2 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.2))
summary = summary(as.mcmc(model.fit.2))
summary$statistics
results2 = as.matrix(as.mcmc(model.fit.2))

## pvalue without claws, without crayfish
p.val(results2[,1] - results2[,2])


data = list(X = with.claw.frags, Y = base, N1 = length(with.claw.frags), N2 = length(without.claw.frags))
model.fit.3 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.3))
summary = summary(as.mcmc(model.fit.3))
summary$statistics
results3 = as.matrix(as.mcmc(model.fit.3))

sink('claws_comparisons.txt', append = TRUE)
cat('###nb.fragments\n\n')

cat('pvalue with claws vs without claws\n')
p.val(results[,1] - results[,2])
cat('\n')

cat('pvalue without claws vs temoin\n')
p.val(results2[,1] - results2[,2])
cat('\n')

cat('pvalue with claws vs temoin\n')
p.val(results3[,1] - results3[,2])
cat('\n')

sink()


##### LDMC ###############
rm(list = setdiff(ls(), lsf.str()))
tab.esp = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/Traits_ESP_LDMC.csv', header = T)
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)
with.claw.LDMC = tab$LDMC[tab$herbi == 'E' & tab$sps == 'Jussie']
without.claw.LDMC = tab.esp$LDMC[tab.esp$sps == 'Jussie']
base = tab$Nbfragments[tab$herbi == 'T' & tab$sps == 'Jussie']

## pinces vs sans pinces ###
data = list(X = with.claw.LDMC, Y = without.claw.LDMC, N1 = length(with.claw.LDMC), N2 = length(without.claw.LDMC))


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

distrib.with.claw = results[,1]
distrib.without.claw = results[,2]

## pvalue with claws, without claws
p.val(results[,1] - results[,2])


data = list(X = without.claw.LDMC, Y = base, N1 = length(without.claw.LDMC), N2 = length(base))
model.fit.2 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.2))
summary = summary(as.mcmc(model.fit.2))
summary$statistics
results2 = as.matrix(as.mcmc(model.fit.2))

## pvalue without claws, without crayfish 
p.val(results2[,1] - results2[,2])


data = list(X = with.claw.LDMC, Y = base, N1 = length(with.claw.LDMC), N2 = length(base))
model.fit.3 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.3))
summary = summary(as.mcmc(model.fit.3))
summary$statistics
results3 = as.matrix(as.mcmc(model.fit.3))


sink('claws_comparisons.txt', append = TRUE)
cat('###LDMC\n\n')

cat('pvalue with claws vs without claws\n')
p.val(results[,1] - results[,2])
cat('\n')

cat('pvalue without claws vs temoin\n')
p.val(results2[,1] - results2[,2])
cat('\n')

cat('pvalue with claws vs temoin\n')
p.val(results3[,1] - results3[,2])
cat('\n')

sink()


########################## Mspi #####################################



sink('claws_comparisons.txt', append = TRUE)
cat('####################\n')
cat('###### Mspi ######\n')
cat('####################\n\n')
sink()
##########333 nb fragments ############
rm(list = setdiff(ls(), lsf.str()))
tab.esp = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/Traits_ESP.csv', header = T, sep = ',')
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv')
with.claw.frags = tab$Nbfragments[tab$herbi == 'E' & tab$sps == 'MSpi']
without.claw.frags = tab.esp$Nbfragments[tab.esp$sps == 'Mspi']
base = tab$Nbfragments[tab$herbi == 'T' & tab$sps == 'MSpi']

## pinces vs sans pinces ###
data = list(X = with.claw.frags, Y = without.claw.frags, N1 = length(with.claw.frags), N2 = length(without.claw.frags))


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

distrib.with.claw = results[,1]
distrib.without.claw = results[,2]

## pvalue with claws, without claws
p.val(results[,1] - results[,2])


data = list(X = without.claw.frags, Y = base, N1 = length(with.claw.frags), N2 = length(without.claw.frags))
model.fit.2 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.2))
summary = summary(as.mcmc(model.fit.2))
summary$statistics
results2 = as.matrix(as.mcmc(model.fit.2))

## pvalue without claws, without crayfish
p.val(results2[,1] - results2[,2])


data = list(X = with.claw.frags, Y = base, N1 = length(with.claw.frags), N2 = length(without.claw.frags))
model.fit.3 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.3))
summary = summary(as.mcmc(model.fit.3))
summary$statistics
results3 = as.matrix(as.mcmc(model.fit.3))

sink('claws_comparisons.txt', append = TRUE)
cat('###nb.fragments\n\n')

cat('pvalue with claws vs without claws\n')
p.val(results[,1] - results[,2])
cat('\n')

cat('pvalue without claws vs temoin\n')
p.val(results2[,1] - results2[,2])
cat('\n')

cat('pvalue with claws vs temoin\n')
p.val(results3[,1] - results3[,2])
cat('\n')


sink()


##### LDMC ###############
rm(list = setdiff(ls(), lsf.str()))
tab.esp = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/Traits_ESP_LDMC.csv', header = T)
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)
with.claw.LDMC = tab$LDMC[tab$herbi == 'E' & tab$sps == 'Mspi']
without.claw.LDMC = tab.esp$LDMC[tab.esp$sps == 'Mspi']
base = tab$Nbfragments[tab$herbi == 'T' & tab$sps == 'Mspi']

## pinces vs sans pinces ###
data = list(X = with.claw.LDMC, Y = without.claw.LDMC, N1 = length(with.claw.LDMC), N2 = length(without.claw.LDMC))


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

distrib.with.claw = results[,1]
distrib.without.claw = results[,2]

## pvalue with claws, without claws
p.val(results[,1] - results[,2])


data = list(X = without.claw.LDMC, Y = base, N1 = length(without.claw.LDMC), N2 = length(base))
model.fit.2 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.2))
summary = summary(as.mcmc(model.fit.2))
summary$statistics
results2 = as.matrix(as.mcmc(model.fit.2))

## pvalue without claws, without crayfish 
p.val(results2[,1] - results2[,2])


data = list(X = with.claw.LDMC, Y = base, N1 = length(with.claw.LDMC), N2 = length(base))
model.fit.3 <- jags.parallel(data = data,
                             model.file = 'Welch_t.test.txt', 
                             # model.file = "repeated_glm.txt"
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 100000,
                             n.burnin = 50000,
                             DIC = FALSE)

plot(as.mcmc(model.fit.3))
summary = summary(as.mcmc(model.fit.3))
summary$statistics
results3 = as.matrix(as.mcmc(model.fit.3))


sink('claws_comparisons.txt', append = TRUE)
cat('###LDMC\n\n')

cat('pvalue with claws vs without claws\n')
p.val(results[,1] - results[,2])
cat('\n')

cat('pvalue without claws vs temoin\n')
p.val(results2[,1] - results2[,2])
cat('\n')

cat('pvalue with claws vs temoin\n')
p.val(results3[,1] - results3[,2])
cat('\n')


sink()
