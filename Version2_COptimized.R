# Must install RapidScore C++ Library First!
#source("~/Documents/Research/SequenceOptimization/SpecificityOptimization.R")
library(RapidScore)
source("~/Documents/Research/SELEX/MultinomialPaper/Submissions/PNAS/Figures/plotting_format.R")
all.models = load.hox.models()

########################
### Define Functions ###
########################
# Convenient barrier function. Height determines strength of penality, sensitivity is how 'square' the 
# logistic response is, and width is how 'wide' the non-asymptotic regions are
barrier.func = function(height, sensitivity, width, cost) {
  # log(1/99) = -4.59512
  return(height/(1+exp(-4.59512*sensitivity/width*cost)) )
}

# Base to number
a.to.n.conversion = c(0:3)
names(a.to.n.conversion) = c("A","C","G","T")
n.to.a.conversion = c("A","C","G","T")

# biased wang-landau 
wang.landau = function(score.object, start.seq, cost.function, verbose, mutation.penalty) {
  # Size of energy divisions
  energy.div.size = 1
  f = 1
  minCount = 100
  minIters = 10000
  fCrit = .9
  tol = .5
  quantile.cutoff=.75
  
  # Convert sequence to base 4 encoding
  seq.len = length(start.seq)
  #seq.pos = 1:seq.len
  seq.pos = 1:(seq.len-1)
  max.mutate.bases = 1 #1:round(seq.len/20)
  xn = rep(0, seq.len)
  for (i in 1:seq.len) {
    xn[i] = a.to.n.conversion[start.seq[i]]
  }
  base.seq = xn
  curr.opt.seq = xn
  
  # Compute starting value
  scores = score.object$scoreBulk(curr.opt.seq)
  nModels = length(scores)
  base.scores = scores
  start.cost = cost.function(scores, base.scores, 0)
  En = start.cost
  curr.opt.cost = En
  if (verbose) {
    cat(paste0("Starting Cost: ", start.cost, "\n"))
  }
  
  if (mutation.penalty!=0) {
    # Compute penalty based on the distribution of (positive) changes 
    score.dist = NULL
    for (i in 1:seq.len) {
      for (curr.base in 0:3) {
        # Skip for the same base
        if (curr.base==base.seq[i]) {
          next
        }
        test.seq = base.seq
        test.seq[i] = curr.base
        scores = score.object$scoreBulk(test.seq)
        test.cost = cost.function(scores, base.scores, 0)-start.cost
        if (test.cost > 0) {
          score.dist = c(score.dist, test.cost)
        }
      }
    }
    # Find percentile value from score distribution and set mutation penalty to it
    mutation.penalty = as.numeric(quantile(score.dist, quantile.cutoff))
    print(mutation.penalty)
  }
  
  # Initialize Wang-Landau scheme
  Emin = floor(min(0, En)/energy.div.size)*energy.div.size
  Emax = ceiling(max(0, En)/energy.div.size)*energy.div.size
  nDivs= ((Emax-Emin)/energy.div.size+1)
  H = rep(0, nDivs)
  S = rep(0, nDivs)
  x0= xn
  E0= En
  I0= round((E0-Emin)/energy.div.size)+1
  iters = 0
  
  # Loop continuously until convergence
  while (TRUE) {
    # Propose a new state
    xn = x0
    mutate.pos = sample(seq.pos, 1)
    xn[mutate.pos:(mutate.pos+1)] = sample(0:3, 2, replace=T)
    # nMuts = sample(max.mutate.bases, 1)
    # mutate.pos = sample(seq.pos, nMuts)
    # xn[mutate.pos] = sample(1:4, nMuts, replace=T)
    # Compute new cost
    scores = score.object$scoreBulk(xn)
    En = cost.function(scores, base.scores, mutation.distance = mutation.penalty*sum(xn!=base.seq))
    # See if the latest sequence is the best one
    if (En > curr.opt.cost) {
      curr.opt.cost = En
      curr.opt.seq = xn
    }
    if (En < Emin) {
      next
    }
    # Compute new index
    In = round((En-Emin)/energy.div.size)+1
    # Update range and restart wang-landau process if a sequence is produced that is out of range
    if (In<=0 || In > nDivs) {
      cat(paste0("Old Range: ", Emin, ":", Emax,"; "))
      # Recompute range
      Emin.new = floor(min(Emin, En)/energy.div.size)*energy.div.size
      Emax.new = ceiling(max(Emax, En)/energy.div.size)*energy.div.size
      cat(paste0("New Range: ", Emin.new, ":", Emax.new))
      # Restart Wang-landau
      if (f==1) {
        # Do not need a hard restart
        minH = min(H[H>0])
        H = c(rep(0, Emin-Emin.new), H, rep(0, Emax.new-Emax))
        S = c(rep(0, Emin-Emin.new), S, rep(0, Emax.new-Emax))
        nDivs= length(H)
        Emin = Emin.new
        Emax = Emax.new
        In = round((En-Emin)/energy.div.size)+1
        cat(paste0("; Total Iterations: ", iters, "; mean: ", mean(H[H>0]), "; sd: ", sd(H[H>0]), "; min: ", minH, "\n"))
      } else {
        # Hard restart required
        Emin = Emin.new
        Emax = Emax.new
        nDivs= ((Emax-Emin)/energy.div.size+1)
        H = rep(0, nDivs)
        S = rep(0, nDivs)
        x0= xn
        E0= En
        I0= round((E0-Emin)/energy.div.size)+1
        iters = 0
        f = 1
        cat("; Restarted sampling!\n")
        next
      }
    }
    # Compute latest acceptance
    a =exp(S[I0] - S[In])
    # Metropolis-Hastings step
    if (runif(1) <= a) {
      # Accept current proposal
      x0 = xn
      E0 = En
      I0 = In
    }
    # Update
    S[I0] = S[I0] + f
    H[I0] = H[I0] + 1
    iters = iters + 1
    # Ensure that all energy bins have been visited an even number of times
    if (min(H[H>0]) > minCount) {
      if (min(H[H>0]) >= fCrit*mean(H[H>0])) {
        cat(paste0("Flatness criteria reached for f = ",f," in ",iters," iterations.\n"))
        f = f/2;
        if (f<=tol || iters>minIters) {
          # Complete. Exit
          break
        }
        iters = 0;
        H = rep(0, nDivs)
      }
    }
  }
  
  # Convert optimal sequence to string representation
  opt.seq = start.seq
  for (i in 1:seq.len) {
    opt.seq[i] = n.to.a.conversion[curr.opt.seq[i]]
  }
  opt.seq = paste0(opt.seq, collapse="")
  return(c(start.cost, opt.seq, curr.opt.cost))
}

# Sequence optimizer. REQUIRES a seed sequence. All scores are automatically adjusted to accomodate the TOTAL sequence score
optimal.sequence = function(..., model.list = NA, cost.function, seq.len = NA, seed.seq,
                            return.max.only = TRUE, verbose = FALSE, penalize.mutations = FALSE) {
  warm.start.multiple = 25
  mutate.dinucs = TRUE
  use.wang.landau = TRUE
  
  if (is.na(model.list)) {
    input.models = list(...)
  } else {
    input.models = model.list
  }
  nModels = length(input.models)
  
  # Create scoring object
  score.object = new(RapidScore)

  # Next, load models
  ks = 0
  for (i in 1:nModels) {
    nuc = input.models[[i]]$NB
    dinuc = input.models[[i]]$DB
    score.object$add(nuc, dinuc)
    ks = max(ks, length(nuc)/4)
  }
  
  # Next, optimize the sequence. Begin by appropriately handling the length of 
  # the string to be optimized and the seed.
  seq.len = max(ks, seq.len, nchar(seed.seq), na.rm = TRUE)
  # Convert seed into a character array
  seed.seq = substring(seed.seq, 1:nchar(seed.seq), 1:nchar(seed.seq))
  # Handle the case where the length of the seed sequence is less than either seq.len or the minimum length requirement [max(ks)]
  if (length(seed.seq)<seq.len) {
    pad.len = (seq.len-length(seed.seq))/2
    # Define length of random field
    l.len = ceiling(pad.len)
    r.len = floor(pad.len)
  } else {
    l.len = 0
    r.len = 0
  }
  
  if (verbose) {
    cat(paste0("Sequence length to optimize: ",seq.len,"; Left Random Field: ", l.len,", Right Random Field: ", r.len, "\n"))
  }
  
  # Set mutation penalty 
  if (!penalize.mutations) {
    mutation.penalty = 0
  } else {
    mutation.penalty = 1
  }
  
  run.outputs = NULL
  # Optimize using wang-landau sampling. Construct starting sequence
  start.seq = c(sample(c("A","C","G","T"), l.len, replace = TRUE), seed.seq, sample(c("A","C","G","T"), r.len, replace = TRUE))
  # Run optimization
  output = wang.landau(score.object, start.seq, cost.function, verbose, mutation.penalty)
  run.outputs = rbind(run.outputs, c(paste(start.seq, collapse=""), output))
  run.outputs = data.frame(run.outputs, stringsAsFactors = FALSE)
  names(run.outputs) = c("StartSequence", "StartCost", "EndSequence", "EndCost")
  
  print(run.outputs)
  
  if (return.max.only) {
    return(run.outputs[which(run.outputs$EndCost==max(run.outputs$EndCost))[1],])
  } else {
    return(run.outputs)
  }
}

#############################################################
### Optimize a variety of sequences in different settings ###
#############################################################
# Load all models
input.models = all.models[[2]][[2]]
model.names = as.character(all.models[[1]]$Protein)
model.names[all.models[[1]]$Dataset=="Dimer"] = paste0("Exd", model.names[all.models[[1]]$Dataset=="Dimer"])
names(input.models) = model.names
for (i in 1:length(input.models)) {
  if (all.models[[1]]$Dataset[i]=="Dimer" || i %in% c(17, 18)) {
    input.models[[i]] = input.models[[i]][[1]]
  }
}

# evaluate all sequences and create plots
genome.profile = function(sequence) {
  info = all.models[[1]]
  for (i in 1:nrow(info)) {
    p = plot.score.genome(genomicSequence = DNAString(sequence), fits = all.models[[2]], index = i, mode = 1)
    print(p + ggtitle(paste0(info$Protein[i], " ", info$Dataset[i])))
  }
}

# simple cost function: Maintain/improve the affinity of the first protein, drop the affinity of the rest
cost.fun = function(scored.list, base.scores, mutation.distance) {
  cost = log(sum(exp(scored.list[[1]]))/sum(exp(base.scores[[1]]))) - 
    barrier.func(1000, 1, exp(max(scored.list[[1]]))/10, 
                 .85*exp(max(base.scores[[1]]))-exp(max(scored.list[[1]])))
  for (i in 2:length(scored.list)) {
    cost = cost + log(sum(exp(base.scores[[i]]))/sum(exp(scored.list[[i]])))
  }
  return(cost-mutation.distance)
}

# Optimize the fkh250 sequence
fkh250 = "CCTCGTCCCACAGCTGGCGATTAATCTTGACATTGAG"
# Loop over all proteins
final.seqs = fkh250
final.costs = 0
final.seqs.p = fkh250
final.costs.p = 0
for (i in 1:length(input.models)) {
  # Select a single model
  idx = 1:length(input.models)
  idx = c(i, idx[-i])
  curr.list = input.models[idx]
  cat(paste0("Optimizing ", names(curr.list)[1], ":\n"))
  out = optimal.sequence(model.list = curr.list, cost.function = cost.fun, seed.seq = fkh250, verbose = TRUE)
  final.seqs = c(final.seqs, out$EndSequence)
  final.costs= c(final.costs,out$EndCost)
  out = optimal.sequence(model.list = curr.list, cost.function = cost.fun, seed.seq = fkh250, verbose = TRUE, penalize.mutations = TRUE)
  final.seqs.p = c(final.seqs.p, out$EndSequence)
  final.costs.p= c(final.costs.p,out$EndCost)
}

# Now compute final affinities 
max.affinity = matrix(data = 0, nrow = length(input.models), ncol=length(final.seqs))
tot.affinity = matrix(data = 0, nrow = length(input.models), ncol=length(final.seqs))
for (i in 1:length(final.seqs)) {
  for (j in 1:length(input.models)) {
    p = score.genome(genomicSequence = DNAString(final.seqs[i]), fits = all.models[[2]], index = j, mode = 1)
    p = p/max.seq(fits = all.models[[2]], index=j, mode=1)$MaxAffinity
    max.affinity[j,i] = max(p)
    tot.affinity[j,i] = sum(p)
  }
}

max.affinity.p = matrix(data = 0, nrow = length(input.models), ncol=length(final.seqs.p))
tot.affinity.p = matrix(data = 0, nrow = length(input.models), ncol=length(final.seqs.p))
for (i in 1:length(final.seqs.p)) {
  for (j in 1:length(input.models)) {
    p = score.genome(genomicSequence = DNAString(final.seqs.p[i]), fits = all.models[[2]], index = j, mode = 1)
    p = p/max.seq(fits = all.models[[2]], index=j, mode=1)$MaxAffinity
    max.affinity.p[j,i] = max(p)
    tot.affinity.p[j,i] = sum(p)
  }
}