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
# 
# pairwise.comps = function(means){
#   cat('---------- Pairwise comparisons: ---------\n')
#   cat('variables \t Effect Size \t p.value\n')
#   for (i in 1:dim(means)[2]){
#     if (i<dim(means)[2]){
#       for (j in (i+1):dim(means)[2]){
#         cat(names(means)[i], ' - ', names(means)[j],':\t', 
#             format(round(mean(means[,i] - means[,j]),3), digits = 3), ' ', 
#             '(', format(round(p.val(means[,i] - means[,j]), digits = 3), 3), ')', '\t',i, ' ',j, '\n')
#       }
#     }
#   }
#   cat('------------------------------------------\n')
# }


pairwise.comps = function(means){
  cat('---------- Pairwise comparisons: ---------\n')
  cat('variables \t Effect Size \t p.value\n')
  for (i in 1:dim(means)[2]){
    if (i<dim(means)[2]){
      for (j in (i+1):dim(means)[2]){
        cat(names(means)[i], ' - ', names(means)[j],':\t', 
            round(mean(means[,i] - means[,j]), 3), ' ', 
            '(', round(p.val(means[,i] - means[,j]), 3), ')', '\t',i, ' ',j, '\n', sep = '')
      }
    }
  }
  cat('------------------------------------------\n')
}


