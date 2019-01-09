rm(list = ls())
library(R2jags)
library(gridExtra)
library(ggplot2)
setwd('/homes/bg33novu/projects/Lise_ecrevisses/traits_and_invasions/')

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


######## definition of plot functions ################


do.data.frame = function(data, order, species, letters){
  means = tapply(data, order, mean, na.rm = TRUE)
  nb = table(!is.na(data), order)
  nb = nb[row.names(nb) == 'TRUE']
  se = abs(tapply(data, order, sd, na.rm = TRUE)/sqrt(nb))
  d = data.frame(y = means, 
                 lower = means - se, 
                 upper= means + se, 
                 pred = pred, 
                 species = species,
                 letters = letters)
  return(d)
}

doplot = function(d, trait, sp.studied, unit, points, preds, species2){
  
  dodge = position_dodge(.3)
  
  ylab = generate.ylab(trait, sp.studied, unit)
  
  fig = ggplot(d, aes(x = pred, y = y, color = species))+
    scale_colour_manual(values = set.colors(species))+
    # geom_boxplot(data = data.frame(preds = preds, points = points, species2 = species2), 
    #              aes(x = preds, y = points, color = species2), 
    #              cex = 0.5, fill = 'grey40', position = dodge)
    geom_errorbar(aes(ymin=lower, ymax=upper, color = species),
                  width = 0.1, size=0.5, position = dodge)+
    geom_point(cex = 2, fill = 'grey40', position = dodge)+
    labs(y = '', x = '')+
    geom_text(aes(x = pred, y = y, label = letters), 
              size = 3,
              nudge_x = c(-0.25, 0.25, -0.25, 0.25),
              colour = 'black',
              show.legend = FALSE)+
    # scale_x_discrete(labels= c('', ''))+
    geom_point(data = data.frame(preds = preds, points = points, species2 = species2), aes(x = preds, y = points, color = species2), cex = 0.5, fill = 'grey40', position = dodge)+
    # geom_boxplot(data = data.frame(preds = preds, points = points, species2 = species2), 
    #                aes(x = preds, y = points, color = species2), 
    #                cex = 0.5, fill = 'grey40', position = dodge)+
    labs(color='Neighbooring species')+
    theme_classic()+
    theme(
      # axis.title.y = element_text(size = 12),
      # axis.title.x = element_text(size = 12),
      # legend.title = element_text(size = 12),
      # legend.text = element_text(size = 10, face = 'italic'),
      axis.text.y = element_text(size = 7),
      legend.position = 'none'
      # axis.text.x = element_blank()
    )
  
  return(fig)
}

generate.ylab = function(trait, sp.studied, unit = NA){
  if (unit == ''){
    return(bquote(paste(.(trait), ' of  ', ~italic(.(sp.studied)), sep = '')))
  }else{
    return(bquote(paste(.(trait), ' of  ', ~italic(.(sp.studied)), ' (', .(unit), ')', sep = '')))
  }
}
set.colors = function(species){
  cols = species
  cols[species == 'L. grandiflora'] = 'red2'
  cols[species == 'M. spicatum'] = 'green4'
  cols[species == 'E. densa'] = 'orange1'
  return(cols)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


pred = factor(c('T', 'T', 'E', 'E'), levels = c('T', 'E'))

########## Jussie ###########
rm(list = setdiff(ls(), lsf.str()))

tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)
Jussie = tab[tab$sps == 'Jussie', ]
attach(Jussie)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)

species = c('E. densa', 'M. spicatum', 'E. densa', 'M. spicatum')
species2 = as.character(Jussie$voisin)
species2[species2 == 'Ms'] = 'M. spicatum'
species2[species2 == 'Eg'] = 'E. densa'
# species2 = factor(species2, order = )
order = factor(interaction(Jussie$ttt,Jussie$voisin), levels = c('T.Eg', 'T.Ms', 'E.Eg', 'E.Ms'))

factors = as.numeric(interaction(voisin,herbi))
tanks = paste('t', ue, sep = '')
tanks = as.numeric(as.factor(as.character(ue)))
nb_tanks = max(tanks)

#### SLA #####
data = list(y = log10(SLA), X = design, N = length(SLA), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.file = "repeated_glm.txt"
cook.values = cook.distance.distrib.b(data, model.file = "repeated_glm.txt")
cook.values[cook.values > 6*mean(cook.values)]

# corresponding plots: 
################# SLA 
jpeg('/homes/bg33novu/projects/Lise_ecrevisses/paper/outliers_tex/jussie_SLA.jpeg')

# layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
dframe = data.frame(cook.values, x = c(1:length(cook.values)))
cooks = ggplot()+
  geom_point(aes(x = dframe$x, y = cook.values)) +
  labs(x = "Index", y = 'Cook values')+
  geom_hline(yintercept = 6*mean(cook.values))+
  geom_hline(yintercept = 4*mean(cook.values), linetype="dashed")+
  theme_classic()

# plot(cook.values, main = 'SLA of L.grandiflora')
# abline(h=6*mean(cook.values))
# abline(h=4*mean(cook.values), lty = 2)

ylimits = c(min(log10(SLA), na.rm = TRUE)*0.9, max(log10(SLA), na.rm = TRUE)*1.1)
unit = ''
trait = 'SLA'
letters = c('ab', 'ab', 'a', 'b')
pred = factor(c('T', 'T', 'E', 'E'), levels = c('T', 'E'))
d = do.data.frame(log10(SLA), order, species, letters)
with = doplot(d, trait, 'L. grandiflora', unit, log10(SLA), ttt, species2) 
with = with  + coord_cartesian(ylim = ylimits) 

letters = c('a', 'b', 'a', 'b')
outliers = which(Jussie$ue == 56)
d.without = do.data.frame(log10(SLA)[-outliers], order[-outliers], species, letters)
without = doplot(d.without, trait, 'L. grandiflora', unit, log10(SLA)[-outliers], ttt[-outliers], species2[-outliers])
without = without + coord_cartesian(ylim = ylimits) 

leg.plot = without + 
  theme(
  legend.position = 'bottom',
  legend.text = element_text(face = 'italic')
)

lay = matrix(c(1,1,2,3), nrow = 2, byrow = TRUE)
legend = g_legend(leg.plot)
grid.arrange(cooks + ggtitle("A)"), 
             with + labs(y = expression(paste('Log10 of Specific',' Leaf area (','mm'^'2','.mg'^'-1',')', sep = '')), x = 'With outliers')
             + ggtitle("B)"), 
             without + labs(y = '', x = 'Without outliers') + ggtitle("C)"), 
             layout_matrix = lay, bottom = legend)
dev.off()

### LDMC ######
data = list(y = log10(LDMC), X = design, N = length(LDMC), nb.factors = ncol(design), id.factor = factors, nb_tanks = nb_tanks, id = tanks)
params = c("beta")
model.file = "repeated_glm.txt"
cook.values = cook.distance.distrib.b(data, model.file = "repeated_glm.txt")
rbind(which(cook.values > 6*mean(cook.values)), cook.values[cook.values > 6*mean(cook.values)])

# corresponding plots: 
plot(cook.values, main = 'LDMC of L.grandiflora')
abline(h=6*mean(cook.values))
abline(h=4*mean(cook.values), lty = 2)


# corresponding plots: 
################# LDMC 
jpeg('/homes/bg33novu/projects/Lise_ecrevisses/paper/outliers_tex/jussie_LDMC.jpeg')

# layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
dframe = data.frame(cook.values, x = c(1:length(cook.values)))
cooks = ggplot()+
  geom_point(aes(x = dframe$x, y = cook.values)) +
  labs(x = "Index", y = 'Cook values')+
  geom_hline(yintercept = 6*mean(cook.values))+
  geom_hline(yintercept = 4*mean(cook.values), linetype="dashed")+
  theme_classic()

ylimits = c(min(log10(LDMC), na.rm = TRUE)*0.9, max(log10(LDMC), na.rm = TRUE)*1.1)
unit = ''
trait = 'LDMC'
letters = c('a', 'a', 'a', 'a')
pred = factor(c('T', 'T', 'E', 'E'), levels = c('T', 'E'))
d = do.data.frame(log10(LDMC), order, species, letters)
with = doplot(d, trait, 'L. grandiflora', unit, log10(LDMC), ttt, species2) 
with = with  + coord_cartesian(ylim = ylimits) 

letters = c('a', 'b', 'ab', 'ab')
outliers = which(Jussie$ue == 56)
outliers = c(outliers, 64)
d.without = do.data.frame(log10(LDMC)[-outliers], order[-outliers], species, letters)
without = doplot(d.without, trait, 'L. grandiflora', unit, log10(LDMC)[-outliers], ttt[-outliers], species2[-outliers])
without = without  + coord_cartesian(ylim = ylimits) 

leg.plot = without + 
  theme(
    legend.position = 'bottom',
    legend.text = element_text(face = 'italic')
  )

lay = matrix(c(1,1,2,3), nrow = 2, byrow = TRUE)
legend = g_legend(leg.plot)
grid.arrange(cooks + ggtitle("A)"), 
             with + labs(y = expression(paste('Log10 of Leaf Dry Matter Content (','mg.g'^'-1',')', sep = '')), x = 'With outliers') + ggtitle("B)"), 
             without + labs(y = '', x = 'Without outliers') + ggtitle("C)"), 
             layout_matrix = lay, bottom = legend)
dev.off()




detach(Jussie)





############# Egeria #########
### LDMC ####
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_leaves.csv', header = T)
Egeria = tab[tab$sps == 'Egeria', ]
attach(Egeria)

species = c('L. grandiflora', 'M. spicatum', 'L. grandiflora', 'M. spicatum')
species2 = as.character(Egeria$voisin)
species2[species2 == 'Ms'] = 'M. spicatum'
species2[species2 == 'Ju'] = 'L. grandiflora'

order = factor(interaction(ttt,voisin), levels = c('T.Ju', 'T.Ms', 'E.Ju', 'E.Ms'))


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
# plot(cook.values, main = 'LDMC of E. densa')
# abline(h=6*mean(cook.values))
# abline(h=4*mean(cook.values), lty = 2)
################# LDMC 
jpeg('/homes/bg33novu/projects/Lise_ecrevisses/paper/outliers_tex/Egeria_LDMC.jpeg')

# layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
dframe = data.frame(cook.values, x = c(1:length(cook.values)))
cooks = ggplot()+
  geom_point(aes(x = dframe$x, y = cook.values)) +
  labs(x = "Index", y = 'Cook values')+
  geom_hline(yintercept = 6*mean(cook.values))+
  geom_hline(yintercept = 4*mean(cook.values), linetype="dashed")+
  theme_classic()

ylimits = c(min(log10(LDMC), na.rm = TRUE)*0.9, max(log10(LDMC), na.rm = TRUE)*1.1)
unit = ''
trait = 'LDMC'
letters = c('a', 'a', 'a', 'a')
pred = factor(c('T', 'T', 'E', 'E'), levels = c('T', 'E'))
d = do.data.frame(log10(LDMC), order, species, letters)
with = doplot(d, trait, 'E. densa', unit, log10(LDMC), ttt, species2) 
with = with  + coord_cartesian(ylim = ylimits) 

letters = c('a', 'ab', 'b', 'ab')
outliers = 16

d.without = do.data.frame(log10(LDMC)[-outliers], order[-outliers], species, letters)
without = doplot(d.without, trait, 'E. densa', unit, log10(LDMC)[-outliers], ttt[-outliers], species2[-outliers])
without = without  + coord_cartesian(ylim = ylimits) 

leg.plot = without + 
  theme(
    legend.position = 'bottom',
    legend.text = element_text(face = 'italic')
  )

lay = matrix(c(1,1,2,3), nrow = 2, byrow = TRUE)
legend = g_legend(leg.plot)
grid.arrange(cooks + ggtitle("A)"), 
             with + labs(y = expression(paste('Log10 of Leaf Dry Matter Content (','mg.g'^'-1',')', sep = '')), x = 'With outliers') + ggtitle("B)"), 
             without + labs(y = '', x = 'Without outliers') + ggtitle("C)"), 
             layout_matrix = lay, bottom = legend)
dev.off()
detach(Egeria)

############# MSpi ################

### RGR #####
rm(list = setdiff(ls(), lsf.str()))
tab = read.csv('/homes/bg33novu/projects/Lise_ecrevisses/data_individuals.csv', header = T)
MSpi = tab[tab$sps == 'MSpi', ]
attach(MSpi)
voisin = droplevels(voisin)
herbi = droplevels(herbi)
design = model.matrix( ~  herbi * voisin)

species = c('E. densa', 'L. grandiflora', 'E. densa', 'L. grandiflora')
species2 = as.character(voisin)
species2[species2 == 'Eg'] = 'E. densa'
species2[species2 == 'Ju'] = 'L. grandiflora'

order = factor(interaction(herbi,voisin), levels = c('T.Eg', 'T.Ju', 'E.Eg', 'E.Ju'))


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
jpeg('/homes/bg33novu/projects/Lise_ecrevisses/paper/outliers_tex/Mspi_RGR.jpeg')

# layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
dframe = data.frame(cook.values, x = c(1:length(cook.values)))
cooks = ggplot()+
  geom_point(aes(x = dframe$x, y = cook.values)) +
  labs(x = "Index", y = 'Cook values')+
  geom_hline(yintercept = 6*mean(cook.values))+
  geom_hline(yintercept = 4*mean(cook.values), linetype="dashed")+
  theme_classic()

ylimits = c(min(RGR, na.rm = TRUE)*1.1, max(RGR, na.rm = TRUE)*1.1)
unit = ''
trait = 'RGR'
letters = c('ab', 'a', 'ab', 'b')
pred = factor(c('T', 'T', 'E', 'E'), levels = c('T', 'E'))
d = do.data.frame(RGR, order, species, letters)
with = doplot(d, trait, 'M. spicatum', unit, RGR, herbi, species2) 
with = with  + coord_cartesian(ylim = ylimits) 

letters = c('ab', 'a', 'ab', 'b')
outliers = 7

d.without = do.data.frame(RGR[-outliers], order[-outliers], species, letters)
without = doplot(d.without, trait, 'M. spicatum', unit, RGR[-outliers], herbi[-outliers], species2[-outliers])
without = without  + coord_cartesian(ylim = ylimits) 

leg.plot = without + 
  theme(
    legend.position = 'bottom',
    legend.text = element_text(face = 'italic')
  )

lay = matrix(c(1,1,2,3), nrow = 2, byrow = TRUE)
legend = g_legend(leg.plot)
grid.arrange(cooks + ggtitle("A)"), 
             with + labs(y = as.expression(bquote('Relative Growth Rate (' ~d^-1~')')), x = 'With outliers') + ggtitle("B)"), 
             without + labs(y = '', x = 'Without outliers') + ggtitle("C)"), 
             layout_matrix = lay, bottom = legend)
dev.off()