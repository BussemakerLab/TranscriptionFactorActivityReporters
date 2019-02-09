# Init for testing
source("~/Documents/Research/SELEX/MultinomialPaper/Submissions/PNAS/Figures/plotting_format.R")
all.models = load.hox.models()
source("~/Documents/Research/SequenceOptimization/SpecificityOptimization.R")

# Optimize a sequence to satisfy an arbitrary model relationship
ExdDfd = all.models[[2]][[2]][[11]][[1]]
ExdLab = all.models[[2]][[2]][[12]][[1]]
ExdScr = all.models[[2]][[2]][[14]][[1]]
ExdUbx = all.models[[2]][[2]][[16]][[1]]

# Define a cost function, send that to sequence optimizer
# Cost function can optimize over scored sequence values (takes in a vector of scored values in ENERGY space)
cost.fun = function(scored.list) {
  # This is a simple one; try to find the sequence that produces the largest difference in affinity
  return(max(scored.list[[1]])-max(scored.list[[2]]))
}

cost.fun2 = function(scored.list) {
  # This is a simple one; try to find the sequence that produces the largest difference in affinity
  return(max(scored.list[[1]])-max(scored.list[[2]]) - barrier.func(20, 1, .05, .05-exp(max(scored.list[[1]]))))
}

# Create aligned and subtracted motifs
mclone = all.models[[2]]
NB = mclone[[2]][[12]][[1]]$NB
DB = mclone[[2]][[12]][[1]]$DB
NB = NB[5:(72-4)]
DB = DB[17:(272-16)]
mclone[[2]][[12]][[1]]$NB = NB
mclone[[2]][[12]][[1]]$DB = DB

NB = mclone[[2]][[14]][[1]]$NB
DB = mclone[[2]][[14]][[1]]$DB
NB = NB[9:72]
DB = DB[33:272]
mclone[[2]][[14]][[1]]$NB = NB
mclone[[2]][[14]][[1]]$DB = DB
ExdLabPrime = mclone[[2]][[12]][[1]]
ExdScrPrime = mclone[[2]][[14]][[1]]

mclone[[2]][[12]][[1]]$NB[1:4] = mclone[[2]][[12]][[1]]$NB[1:4]-as.numeric(log(max.seq(mclone, 12, 1)[1]))
mclone[[2]][[14]][[1]]$NB[1:4] = mclone[[2]][[14]][[1]]$NB[1:4]-as.numeric(log(max.seq(mclone, 14, 1)[1]))
ExdLabExdScr = mclone[[2]][[14]][[1]]
ExdLabExdScr$NB = mclone[[2]][[12]][[1]]$NB - mclone[[2]][[14]][[1]]$NB
ExdLabExdScr$DB = mclone[[2]][[12]][[1]]$DB - mclone[[2]][[14]][[1]]$DB
mclone[[2]][[15]][[1]] = ExdLabExdScr

sghelper = function(seq, idx1, idx2) {
  return(rbind(log(score.genome(genomicSequence = DNAString(seq), fits = mclone, index = idx1, mode = 1)), 
               log(score.genome(genomicSequence = DNAString(seq), fits = mclone, index = idx2, mode = 1))))
}