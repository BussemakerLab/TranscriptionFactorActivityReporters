source("~/Documents/Research/SELEX/MultinomialPaper/Submissions/PNAS/Figures/plotting_format.R")
all.models = load.hox.models()
source("~/Documents/Research/SequenceOptimization/OldVersions/SpecificityOptimization.R")

vvl.seq = "GGGATTTGATTAATTTATTACCG"

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

performance = function(seq) {
  cat("Optimized Lab Affinity: ", final.affinity(seq, 12), "; Fold Change: ", final.affinity(seq, 12)/exp(lab.orig), "\n")
  cat("Optimized Dfd Affinity: ", final.affinity(seq, 11), "; Fold Change: ", final.affinity(seq, 11)/exp(dfd.orig), "\n")
  cat("Optimized Scr Affinity: ", final.affinity(seq, 14), "; Fold Change: ", final.affinity(seq, 14)/exp(scr.orig), "\n")
  cat("Optimized Ubx Affinity: ", final.affinity(seq, 16), "; Fold Change: ", final.affinity(seq, 16)/exp(ubx.orig), "\n")
}

lab.orig = init.energy(12)
dfd.orig = init.energy(11)
scr.orig = init.energy(14)
ubx.orig = init.energy(16)

cost.fun.lab = function(scored.list, null.input) {
  # Find a sequence that improves the binding of labial, reduces the binding of dfd, scr, and ubx, and ensure that labial affinity does not drop below lab.orig
  return((max(scored.list[[1]])-lab.orig) + 
           max((dfd.orig-max(scored.list[[2]])), floor) + 
           max((scr.orig-max(scored.list[[3]])), floor) + 
           max((ubx.orig-max(scored.list[[4]])), floor) - 
           barrier.func(40, 10, exp(lab.orig), exp(lab.orig)-exp(max(scored.list[[1]]))))
}

cost.fun.dfd = function(scored.list, null.input) {
  # Find a sequence that improves the binding of dfd, reduces the binding of lab, scr, and ubx, and ensure that dfd affinity does not drop below lab.orig
  return((max(scored.list[[2]])-dfd.orig)*17 + 
           (lab.orig-max(max(scored.list[[1]]), floor)) + 
           (scr.orig-max(max(scored.list[[3]]), floor)) + 
           (ubx.orig-max(max(scored.list[[4]]), floor)) - 
           barrier.func(40, 10, exp(dfd.orig), exp(dfd.orig)-exp(max(scored.list[[2]]))) -
           barrier.func(10, 10, exp(lab.orig), -exp(lab.orig)+exp(max(scored.list[[1]]))) - 
           barrier.func(10, 10, exp(scr.orig), -exp(dfd.orig)+exp(max(scored.list[[3]]))) - 
           barrier.func(10, 10, exp(ubx.orig), -exp(ubx.orig)+exp(max(scored.list[[4]]))))
}

cost.fun.scr = function(scored.list, null.input) {
  # Find a sequence that improves the binding of scr, reduces the binding of lab, dfd, and ubx, and ensure that scr affinity does not drop below scr.orig
  return((max(scored.list[[3]])-scr.orig)*15 + 
           (lab.orig-max(max(scored.list[[1]]), floor)) + 
           (dfd.orig-max(max(scored.list[[2]]), floor)) + 
           (ubx.orig-max(max(scored.list[[4]]), floor)) - 
           barrier.func(40, 10, exp(scr.orig), exp(scr.orig)-exp(max(scored.list[[3]]))) - 
           barrier.func(10, 10, exp(lab.orig), -exp(lab.orig)+exp(max(scored.list[[1]]))) - 
           barrier.func(10, 10, exp(dfd.orig), -exp(dfd.orig)+exp(max(scored.list[[2]]))) - 
           barrier.func(10, 10, exp(ubx.orig), -exp(ubx.orig)+exp(max(scored.list[[4]]))))
}

cost.fun.ubx = function(scored.list, null.input) {
  # Find a sequence that improves the binding of ubx, reduces the binding of lab, dfd, and scr, and ensure that scr affinity does not drop below ubx.orig
  return((max(scored.list[[4]])-ubx.orig)*5 + 
           (lab.orig-max(max(scored.list[[1]]), floor)) + 
           (dfd.orig-max(max(scored.list[[2]]), floor)) + 
           (scr.orig-max(max(scored.list[[3]]), floor)) - 
           barrier.func(40, 10, exp(ubx.orig), exp(ubx.orig)-exp(max(scored.list[[4]]))))
}

# Define floor
floor = -9
optimal.sequence(ExdLab, ExdDfd, ExdScr, ExdUbx, seed.seq = vvl.seq, cost.function = cost.fun.lab)
performance("GAAAATTGATCGATCGGGACGCG")

optimal.sequence(ExdLab, ExdDfd, ExdScr, ExdUbx, seed.seq = vvl.seq, cost.function = cost.fun.dfd)
performance("GGCGGAATATTAATGATTTTTAG")

optimal.sequence(ExdLab, ExdDfd, ExdScr, ExdUbx, seed.seq = vvl.seq, cost.function = cost.fun.scr)
performance("GACGAATGATAAATTGCTGCGAG")

optimal.sequence(ExdLab, ExdDfd, ExdScr, ExdUbx, seed.seq = vvl.seq, cost.function = cost.fun.ubx)
performance(seq = "GACGGATAATTTACGACCCCCAG")
