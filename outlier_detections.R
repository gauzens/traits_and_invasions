rm(list = ls())
library(R2jags)

##### definition of function computing cook distances ###########

cooks.distance.point.b = function(data, predict, predict.witouts, to.remove, ...){
  MSE = 0
  MSE = MSE + sum((data - predict)^2, na.rm = TRUE)
  MSE = MSE / length(data[!is.na(data)])
  cook = 0
  preds = predict[-to.remove]
  cook = sum((preds - predict.witouts)^2, na.rm = TRUE)
  cook = cook/(MSE * 4)
  return(cook)
}

cook.distance.distrib.b = function(data, ...){
  params = c("beta")
  main.model = jags.parallel(data = data,
                             # model.file = 'repeated_anova.txt'
                             model.file = model.file,
                             parameters.to.save = params,
                             n.chains = 4,
                             n.iter = 70000,
                             n.burnin = 50000,
                             DIC = FALSE)
  results = as.matrix(as.mcmc(main.model))
  results = colMeans(results)
  predicts = data$X %*% results
  
  
  cooks.distrib = rep(NA, length(data$y))
  for (i in 1:length(data$y)){
    cat(i, ' ')
    flush.console()
    data2 = data
    data2$y = data$y[-i]
    data2$X = data$X[-i,] 
    data2$N = data$N - 1
    data2$id.factor = data$id.factor[-i]
    data2$id = data$id[-i]
    # print(data2)
    params = c("beta")
    model.without = jags.parallel(data = data2,
                                  # model.file = 'repeated_anova.txt'
                                  model.file = "repeated_glm.txt",
                                  parameters.to.save = params,
                                  n.chains = 4,
                                  n.iter = 70000,
                                  n.burnin = 50000,
                                  DIC = FALSE)
    results = as.matrix(as.mcmc(model.without))
    results = colMeans(results)
    predicts.without = data2$X %*% results
    cooks.distrib[i] = cooks.distance.point.b(data$y, predicts, predicts.without, i)
  }
  return(cooks.distrib)
}






########## Jussie ###########
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('data_leaves.csv', header = T)
Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)

#### SLA
data = list(y = log10(SLA), X = design, N = length(SLA), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.file = "repeated_glm.txt"
cook.values = cook.distance.distrib.b(data, model.file = "repeated_glm.txt")
cook.values[cook.values > 6*mean(cook.values)]

# corresponding plots: 
plot(cook.values, main = 'SLA of L.grandiflora')
abline(h=6*mean(cook.values))
abline(h=4*mean(cook.values), lty = 2)


### LDMC
data = list(y = log10(LDMC), X = design, N = length(LDMC), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.file = "repeated_glm.txt"
cook.values = cook.distance.distrib.b(data, model.file = "repeated_glm.txt")
rbind(which(cook.values > 6*mean(cook.values)), cook.values[cook.values > 6*mean(cook.values)])

# corresponding plots: 
plot(cook.values, main = 'LDMC of L.grandiflora')
abline(h=6*mean(cook.values))
abline(h=4*mean(cook.values), lty = 2)
detach(Jussie)


############# Egeria #########
### LDMC
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('data_leaves.csv', header = T)
Egeria = tab[tab$sps == 'Egeria', ]
attach(Egeria)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(LDMC), X = design, N = length(LDMC), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.file = "repeated_glm.txt"
cook.values = cook.distance.distrib.b(data, model.file = "repeated_glm.txt")
rbind(which(cook.values > 6*mean(cook.values)), cook.values[cook.values > 6*mean(cook.values)])

# corresponding plots: 
plot(cook.values, main = 'LDMC of E. densa')
abline(h=6*mean(cook.values))
abline(h=4*mean(cook.values), lty = 2)

detach(Egeria)

############# MSpi
### LDMC
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('data_leaves.csv', header = T)
MSpi = tab[tab$sps == 'Mspi', ]
attach(MSpi)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = log10(LDMC), X = design, N = length(LDMC), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.file = "repeated_glm.txt"
cook.values = cook.distance.distrib.b(data, model.file = "repeated_glm.txt")
rbind(which(cook.values > 6*mean(cook.values)), cook.values[cook.values > 6*mean(cook.values)])

# corresponding plots: 
plot(cook.values, main = 'LDMC of M. spicatum')
abline(h=6*mean(cook.values))
abline(h=4*mean(cook.values), lty = 2)
detach(MSpi)

### RGR
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('data_individuals.csv', header = T)
MSpi = tab[tab$sps == 'MSpi', ]
attach(MSpi)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)
data = list(y = RGR, X = design, N = length(RGR), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.file = "fixed_glm.txt"
cook.values = cook.distance.distrib.b(data, model.file = "fixed_glm.txt")
rbind(which(cook.values > 6*mean(cook.values)), cook.values[cook.values > 6*mean(cook.values)])

# corresponding plots: 
plot(cook.values, main =' RGR of M. spicatum')
abline(h=6*mean(cook.values))
abline(h=4*mean(cook.values), lty = 2)
detach(MSpi)
