source("~/Documents/Research/SELEX/MultinomialPaper/Submissions/PNAS/Figures/plotting_format.R")
all.models = load.hox.models()
source("~/Documents/Research/SequenceOptimization/SpecificityOptimization.R")

vvl.seq = "GGGATTTGATTAATTTATTA"

# Optimize a sequence to satisfy an arbitrary model relationship
ExdDfd = all.models[[2]][[2]][[11]][[1]]
ExdLab = all.models[[2]][[2]][[12]][[1]]
ExdScr = all.models[[2]][[2]][[14]][[1]]
ExdUbx = all.models[[2]][[2]][[16]][[1]]

init.energy = function(index) {
  return(log(max(score.genome(genomicSequence = DNAString(vvl.seq), fits = all.models[[2]], index, 1)/as.numeric(max.seq(all.models[[2]], index, 1)[1]))))
}

final.affinity = function(seq, index) {
  return(max(score.genome(genomicSequence = DNAString(seq), fits = all.models[[2]], index, 1)/as.numeric(max.seq(all.models[[2]], index, 1)[1])))
}

lab.orig = init.energy(12)
dfd.orig = init.energy(11)
scr.orig = init.energy(14)
ubx.orig = init.energy(16)

cost.fun3 = function(scored.list) {
  # Find a sequence that improves the binding of labial, reduces the binding of dfd, scr, and ubx, and ensure that labial affinity does not drop below lab.orig
  return((max(scored.list[[1]])-lab.orig) + 
           (dfd.orig-max(scored.list[[2]])) + 
           (scr.orig-max(scored.list[[3]])) + 
           (ubx.orig-max(scored.list[[4]])) - 
           barrier.func(20, 1, exp(lab.orig), exp(lab.orig)-exp(max(scored.list[[1]]))))
}

cost.fun4 = function(scored.list) {
  # Find a sequence that improves the binding of labial, reduces the binding of dfd, scr, and ubx, and ensure that labial affinity does not drop below lab.orig
  return((max(scored.list[[1]])-dfd.orig) + 
           (lab.orig-max(scored.list[[2]])) + 
           (scr.orig-max(scored.list[[3]])) + 
           (ubx.orig-max(scored.list[[4]])) - 
           barrier.func(20, 1, exp(dfd.orig), exp(dfd.orig)-exp(max(scored.list[[1]]))))
}

optimal.sequence(ExdLab, ExdDfd, ExdScr, ExdUbx, seed.seq = vvl.seq, cost.function = cost.fun3)
cat("Optimized Lab Affinity: ", final.affinity("CCCACATGATGGATGGGAGA", 12), "; Fold Change: ", final.affinity("CCCACATGATGGATGGGAGA", 12)/exp(lab.orig), "\n")
cat("Optimized Dfd Affinity: ", final.affinity("CCCACATGATGGATGGGAGA", 11), "; Fold Change: ", final.affinity("CCCACATGATGGATGGGAGA", 11)/exp(dfd.orig), "\n")
cat("Optimized Scr Affinity: ", final.affinity("CCCACATGATGGATGGGAGA", 14), "; Fold Change: ", final.affinity("CCCACATGATGGATGGGAGA", 14)/exp(scr.orig), "\n")
cat("Optimized Ubx Affinity: ", final.affinity("CCCACATGATGGATGGGAGA", 16), "; Fold Change: ", final.affinity("CCCACATGATGGATGGGAGA", 16)/exp(ubx.orig), "\n")

cat("Optimized Lab Affinity: ", final.affinity("ACAAAATGATGGATGGGTTA", 12), "; Fold Change: ", final.affinity("ACAAAATGATGGATGGGTTA", 12)/exp(lab.orig), "\n")
cat("Optimized Dfd Affinity: ", final.affinity("ACAAAATGATGGATGGGTTA", 11), "; Fold Change: ", final.affinity("ACAAAATGATGGATGGGTTA", 11)/exp(dfd.orig), "\n")
cat("Optimized Scr Affinity: ", final.affinity("ACAAAATGATGGATGGGTTA", 14), "; Fold Change: ", final.affinity("ACAAAATGATGGATGGGTTA", 14)/exp(scr.orig), "\n")
cat("Optimized Ubx Affinity: ", final.affinity("ACAAAATGATGGATGGGTTA", 16), "; Fold Change: ", final.affinity("ACAAAATGATGGATGGGTTA", 16)/exp(ubx.orig), "\n")
