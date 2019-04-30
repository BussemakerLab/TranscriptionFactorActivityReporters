# Load depending on OS
if (as.character(Sys.info()['sysname'])=="Windows") {
  # Sequence Engineering Library
  source("~/Columbia Work/Research/SequenceOptimization/base/WLDOS_SequenceOptimizer.R")
  # Load Hox Models
  source("~/Columbia Work/Research/SELEX/MultinomialPaper/Submissions/PNAS/Figures/plotting_format.R")
} else {
  # Sequence Engineering Library
  source("~/Documents/Research/SequenceOptimization/base/WLDOS_SequenceOptimizer.R")
  # Load Hox Models
  source("~/Documents/Research/SELEX/MultinomialPaper/Submissions/PNAS/Figures/plotting_format.R")
}
all.models = load.hox.models()

# Remove UbxIa
input.models = all.models[[2]][[2]]
ubx.idx = grep("UbxIa",all.models[[1]]$Protein)
all.models[[1]] = all.models[[1]][-ubx.idx,]
input.models = input.models[-ubx.idx]
all.models[[2]][[2]] = input.models
all.models[[2]][[1]] = all.models[[2]][[1]][-ubx.idx,]
model.names = as.character(all.models[[1]]$Protein)
model.names[all.models[[1]]$Dataset=="Dimer"] = paste0("Exd", model.names[all.models[[1]]$Dataset=="Dimer"])
names(input.models) = model.names

# Need to adjust the energy scale of all models so that the energy of the maximum sequence is 0. This allows setting a 'nonspecific binding' limit
for (i in 1:length(input.models)) {
  if (!is.na(all.models[[1]]$ModeIdx[i])) {
    input.models[[i]] = input.models[[i]][[1]]
  }
  energy.adjustment = log(max.seq(all.models[[2]], i, 1)$MaxAffinity)
  input.models[[i]]$NB[1:4] = input.models[[i]]$NB[1:4] - energy.adjustment
}

#######################################
### Define Different Cost Functions ###
#######################################
# Simple cost function: Maintain/improve the affinity of the first protein, drop the affinity of the rest
cost.fun = function(scored.list, base.scores) {
  cost = log(sum(exp(scored.list[[1]]))/sum(exp(base.scores[[1]])))
  for (i in 2:length(scored.list)) {
    cost = cost + log(sum(exp(base.scores[[i]]))/sum(exp(scored.list[[i]])))
  }
  return(cost)
}

# Max distance + bonus for significant improvement
max.dist.bonus = function(scored.list, base.scores) {
  top.aff.prot1 = exp(max(scored.list[[1]])-max(base.scores[[1]]))
  
  # Important to keep the affinity of the first protein the same
  if (top.aff.prot1>=1) {
    cost = log10(top.aff.prot1)*10
  } else {
    cost = (log10(top.aff.prot1) - 1)*100
  }
  sum.score.ratio = Inf
  max.score.ratio = Inf
  for (i in 2:length(scored.list)) {
    sum.ratio = sum(exp(base.scores[[i]]))/sum(exp(scored.list[[i]]))
    max.ratio = exp(max(base.scores[[i]])-max(scored.list[[i]]))
    if (sum.ratio>=10) {
      cost = cost + 1
    }
    if (max.ratio>=10) {
      cost = cost + 1
    }
    sum.score.ratio = min(sum.score.ratio, sum.ratio)
    max.score.ratio = min(max.score.ratio, max.ratio)
  }
  if (sum.score.ratio>=1) {
    cost = cost + log10(sum.score.ratio)*10
  } else {
    cost = cost + (log10(sum.score.ratio) - 1)*100
  }
  if (max.score.ratio>=1) {
    cost = cost + log10(max.score.ratio)*10
  } else {
    cost = cost + (log10(max.score.ratio) - 1)*100
  }
  # 'Regularize' the negative cost space. We obviously dont want these solutions, 
  # but it will compress the output and still allow some degree of exploration
  if (cost<0) {
    cost = -log2(-cost)*2
  }
  return(cost)
}

max.improvement = function(scored.list, base.scores) {
  # Improve overall binding of first protein
  prot1.ratio = log10(sum(exp(scored.list[[1]]))/sum(exp(base.scores[[1]])))*100
  if (prot1.ratio>=0) {
    cost = prot1.ratio
  } else {
    cost = prot1.ratio - 100
  }
  
  # Drop binding of the others
  for (i in 2:length(scored.list)) {
    protN.ratio = log10(sum(exp(base.scores[[i]]))/sum(exp(scored.list[[i]])))
    if (protN.ratio>=0) {
      cost = cost + protN.ratio
    } else {
      cost = cost + protN.ratio*10 - 10
    }
  }
  
  # 'Regularize' the negative cost space. We obviously dont want these solutions, 
  # but it will compress the output and still allow some degree of exploration
  if (cost<0) {
    cost = -log2(-cost)*4
  }
  return(cost)
}

#########################
### Run Optimizations ###
#########################
# Optimize the fkh250 sequence
fkh250 = "CCTCGTCCCACAGCTGGCGATTAATCTTGACATTGAG"
# Select a single model
i=14
idx = 1:length(input.models)
idx = c(i, idx[-i])
curr.list = input.models[idx]
cat(paste0("Optimizing ", names(curr.list)[1], ":\n"))
max.improvement.out = optimal.sequence(model.list = curr.list, cost.function = max.improvement, seed.seq = fkh250, verbose = TRUE)