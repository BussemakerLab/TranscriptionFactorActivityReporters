# Must install RapidScore C++ Library First!
library(RapidScore)
library(abind)
source("~/Documents/Research/SELEX/MultinomialPaper/Submissions/PNAS/Figures/plotting_format.R")
all.models = load.hox.models()

########################
### Define Functions ###
########################
# Base to number
a.to.n.conversion = c(0:3)
names(a.to.n.conversion) = c("A","C","G","T")
n.to.a.conversion = c("A","C","G","T")

wang.landau = function(score.object, start.seq, cost.function, verbose=TRUE) {
  # Size of energy divisions and other parameters - play with these?
  energy.div.size = 1
  f = 1
  minCount = 100
  minIters = 1E4
  minLoops = 5
  tol = .8
  globalIterations = 0
  resetIterations = 250000
  
  # Convert sequence to base 4 encoding, and define positions to mutate
  seq.len = length(start.seq)
  seq.pos = 1:seq.len
  xn = rep(0, seq.len)
  for (i in 1:seq.len) {
    xn[i] = a.to.n.conversion[start.seq[i]]
  }
  base.seq = xn
  x0 = xn
  
  # Compute starting values
  scores = score.object$scoreBulk(base.seq)
  base.scores = scores
  start.cost = cost.function(scores, base.scores)
  En = start.cost
  E0 = En
  if (verbose) {
    cat(paste0("Starting Cost: ", start.cost, "\n"))
  }
  
  # Initialize Wang-Landau scheme; start by defining ranges
  Emin = floor(min(0, En)/energy.div.size)*energy.div.size
  Emax = ceiling(max(0, En)/energy.div.size)*energy.div.size
  MutMin = 0
  MutMax = 2
  # Define bin size
  eDivs= ((Emax-Emin)/energy.div.size+1)
  mDivs = MutMax+1
  # State variables; HIdx is a binary matrix storing the indicies of visited bins
  H = matrix(data = 0, ncol = eDivs, nrow = mDivs)
  S = H
  HIdx = H>0
  # Last visited seq array
  last.seq = array(dim = c(mDivs, eDivs, seq.len))
  # List of optimal sequences by mutational distance
  opt.seqs = data.frame(Mutations=MutMin:MutMax, OptimalCost=-Inf)
  opt.seqs$OptimalSeq=list(NA)
  opt.seqs$OptimalCost[1] = E0
  opt.seqs$OptimalSeq[[1]]= x0
  # Monitor previous cycles
  prevOptCost = NA
  prevHIdx = HIdx
  prevNActiveBins = 0
  stableCost = FALSE
  
  # Initialize Wang-Landau Scheme
  IE0= round((E0-Emin)/energy.div.size)+1
  IM0= 1  # Initial bin for edit distance is 0 
  S[IM0, IE0] = S[IM0, IE0] + f
  H[IM0, IE0] = H[IM0, IE0] + 1
  HIdx[IM0, IE0] = TRUE
  prevHIdx[IM0, IE0] = TRUE
  last.seq[IM0, IE0, ] = x0
  iters = 1
  loops = 0
  accept = 0
  reject = 0
  
  # Loop continuously until convergence
  while (TRUE) {
    # Propose a new state with a nucleotide mutation. On special steps, move to a 
    # bin with fewer visits (hopefully). This biased move needs to be compensated
    if (runif(1) <= .005) {   # TODO: does this rate need to scale with sequence length?
      # Compute bin sampling bias distribution
      Hp = H
      Hp[HIdx] = Hp[HIdx]-min(Hp[HIdx])
      Hp[prevHIdx] = Hp[prevHIdx] + 1
      cts = c(Hp)
      cts[cts>0] = sum(cts)/cts[cts>0]
      cts = cts/sum(cts)
      idx = arrayInd(sample(length(H), 1, prob = cts), dim(H))
      xn = last.seq[idx[1],idx[2],]
      xn[sample(seq.pos,1)] = sample(0:3,1)
      # This transition ratio is to account for unbalanced transition probabilities
      trans.ratio = cts[(IE0-1)*mDivs+IM0]/cts[(idx[2]-1)*mDivs+idx[1]]
    } else {
      # 'Standard' move where the sequence is mutated at a single position
      xn = x0
      xn[sample(seq.pos, 1)] = sample(0:3, 1)
      trans.ratio = 1
    }
    
    # Compute new cost
    scores = score.object$scoreBulk(xn)
    En = cost.function(scores, base.scores)
    # Compute mutation distance; Reverse complement symmetry must be accounted for when computing
    Mn = min(sum(xn!=base.seq), sum(rev(-1*(xn-3))!=base.seq)) 
    # Compute new cost and mutation indexes
    IEn = round((En-Emin)/energy.div.size)+1
    IMn = Mn+1
    globalIterations = globalIterations + 1
    
    # Begin by padding the H and S matrices with mutation range increases
    if (IMn > mDivs) {
      # Update and add rows to opt.seq
      for (k in (MutMax+1):(IMn-1)) {
        opt.seqs = rbind(opt.seqs, rep(c(k, -Inf, list(NA))))
      }
      # Recompute Range and divs
      MutMax = IMn-1
      mDivs = IMn
      last.seq = abind(last.seq, array(dim=c(mDivs-nrow(H), ncol(H), seq.len)), along = 1)
      H = rbind(H, matrix(data = 0, ncol = ncol(H), nrow = mDivs-nrow(H)))
      S = rbind(S, matrix(data = 0, ncol = ncol(S), nrow = mDivs-nrow(S)))
      HIdx = rbind(HIdx, matrix(data = FALSE, ncol = ncol(HIdx), nrow = mDivs-nrow(HIdx)))
      prevHIdx = rbind(prevHIdx, matrix(data = FALSE, ncol = ncol(prevHIdx), nrow = mDivs-nrow(prevHIdx)))
    }
    # Check cost ranges; pad H and S with energy range increases
    if (IEn<=0 || IEn > eDivs) {
      # Recompute range and divs
      Emin.new = floor(min(Emin, En)/energy.div.size)*energy.div.size
      Emax.new = ceiling(max(Emax, En)/energy.div.size)*energy.div.size
      H = cbind(matrix(data = 0, ncol = Emin-Emin.new, nrow = nrow(H)), H, matrix(data = 0, ncol = Emax.new-Emax, nrow = nrow(H)))
      S = cbind(matrix(data = 0, ncol = Emin-Emin.new, nrow = nrow(S)), S, matrix(data = 0, ncol = Emax.new-Emax, nrow = nrow(S)))
      HIdx = cbind(matrix(data = FALSE, ncol = Emin-Emin.new, nrow = nrow(HIdx)), HIdx, matrix(data = FALSE, ncol = Emax.new-Emax, nrow = nrow(HIdx)))
      prevHIdx = cbind(matrix(data = FALSE, ncol = Emin-Emin.new, nrow = nrow(prevHIdx)), 
                       prevHIdx, matrix(data = FALSE, ncol = Emax.new-Emax, nrow = nrow(prevHIdx)))
      last.seq = abind(array(dim=c(nrow(H), Emin-Emin.new, seq.len)), last.seq, array(dim=c(nrow(H), Emax.new-Emax, seq.len)), along=2)
      eDivs = ncol(H)
      Emin = Emin.new
      Emax = Emax.new
      IEn = round((En-Emin)/energy.div.size)+1
    }
    
    # See if the latest sequence is the best one; if so, update list of sequences
    if (En > opt.seqs$OptimalCost[IMn]) {
      opt.seqs$OptimalCost[IMn] = En
      opt.seqs$OptimalSeq[[IMn]] = xn
    }
    # Compute latest acceptance
    a =exp(S[IM0, IE0] - S[IMn, IEn])*trans.ratio
    # Metropolis-Hastings step
    if (runif(1) <= a) {
      # Accept current proposal
      x0 = xn
      E0 = En
      IE0 = IEn
      IM0 = IMn
      # Store last accepted sequence
      last.seq[IM0, IE0, ] = x0
      HIdx[IM0, IE0] = TRUE
      prevHIdx[IM0, IE0] = TRUE
      accept = accept + 1
    } else {
      reject = reject + 1
    }
    # Update
    S[IM0, IE0] = S[IM0, IE0] + f
    H[IM0, IE0] = H[IM0, IE0] + 1
    iters = iters + 1
    
    # Monitoring
    if (verbose && iters %% 250000 ==0) {
      cat(paste0("Total Iterations: ", iters, "; mean: $", mean(H[prevHIdx]), "; sd: ", sd(H[prevHIdx]), "; min: ^", 
                 min(H[prevHIdx]), "; active bins: &", sum(prevHIdx), "; empty bins: *", sum(H[prevHIdx]==0), 
                 "; Emin-Emax: ", Emin,"-",Emax,"; accept-reject ratio: @",accept, "\n"))
      accept = 0
      reject = 0
    }
    
    if(iters %% resetIterations == 0 && sum(prevHIdx[!HIdx]==1)>0) {
      S = S*0
      H = H*0
    }
  
    # Ensure that all energy bins have been visited an even number of times
    if (min(H[prevHIdx]) > minCount && min(H[prevHIdx]) >= tol*mean(H[prevHIdx]) && iters>minIters) {
      loops = loops + 1;
      if (verbose) {
        cat(paste0("Flatness criteria reached in ",iters," iterations. ", sum(prevHIdx), " active bins; ", loops, " loops completed.\n"))
      }
      # Terminate if the minimum number of loops is met and if the optimal cost vector has remained constant for 3 loops
      print(prevOptCost)
      print(opt.seqs$OptimalCost)
      if (loops>=minLoops && sum(prevOptCost!=opt.seqs$OptimalCost)==0) {
        cat(paste0("Finished in ", globalIterations, " total iterations.\n"))
        break
      }
      #   cat("Cost is stable.\n")
      #   print(stableCost)
      #   print(loops>=minLoops && stableCost)
      #   if (loops>=minLoops && stableCost) {
      #     cat(paste0("Finished in ", globalIterations, " total iterations.\n"))
      #     break  
      #   } else {
      #     stableCost = TRUE
      #   }
      # } else {
      #   stableCost = FALSE
      # }
      # print(stableCost)
      # Store this round's optimal cost and active bin indices
      prevOptCost = opt.seqs$OptimalCost
      # Reset the sampling for the next round
      H = matrix(data = 0, ncol = eDivs, nrow = mDivs)
      S = H
      HIdx = H>0
      #last.seq = array(dim = c(mDivs, eDivs, seq.len))
      last.seq[IM0, IE0, ] = x0
      S[IM0, IE0] = S[IM0, IE0] + f
      H[IM0, IE0] = H[IM0, IE0] + 1
      HIdx[IM0, IE0] = TRUE
      prevHIdx[IM0, IE0] = TRUE
      globalIterations = globalIterations + iters;
      iters = 1;
    }
  }
  
  # Convert optimal sequence to string representation
  for (i in 1:nrow(opt.seqs)) {
    curr.seq = opt.seqs$OptimalSeq[[i]]
    if (sum(is.na(curr.seq))>0) {
      next
    }
    out.seq = start.seq
    for (j in 1:seq.len) {
      out.seq[j] = n.to.a.conversion[curr.seq[j]+1]
    }
    out.seq = paste0(out.seq, collapse="")
    opt.seqs$OptimalSeq[[i]] = out.seq
  }
  opt.seqs$OptimalSeq = unlist(opt.seqs$OptimalSeq)
  print(opt.seqs)
  return(opt.seqs)
}

# Sequence optimizer. REQUIRES a seed sequence. All scores are automatically adjusted to accomodate the TOTAL sequence score
optimal.sequence = function(..., model.list = NA, cost.function, seq.len = NA, seed.seq, affinity.cutoff = 1E-5, verbose = FALSE) {
  if (is.na(model.list)) {
    input.models = list(...)
  } else {
    input.models = model.list
  }
  nModels = length(input.models)
  
  # Create scoring object
  score.object = new(RapidScore, log(affinity.cutoff))
  
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
  
  # Optimize using wang-landau sampling. Construct starting sequence
  start.seq = c(sample(c("A","C","G","T"), l.len, replace = TRUE), seed.seq, sample(c("A","C","G","T"), r.len, replace = TRUE))
  # Run optimization
  output = wang.landau(score.object, start.seq, cost.function, verbose = verbose)
  return(output)
}

#############################################################
### Optimize a variety of sequences in different settings ###
#############################################################
# Load all models
input.models = all.models[[2]][[2]]
model.names = as.character(all.models[[1]]$Protein)
model.names[all.models[[1]]$Dataset=="Dimer"] = paste0("Exd", model.names[all.models[[1]]$Dataset=="Dimer"])
names(input.models) = model.names
# Need to adjust the energy scale of all models so that the energy of the maximum sequence is 0. This allows setting a 'nonspecific binding' limit
for (i in 1:length(input.models)) {
  if (all.models[[1]]$Dataset[i]=="Dimer" || i %in% c(17, 18)) {
    input.models[[i]] = input.models[[i]][[1]]
  }
  energy.adjustment = log(max.seq(all.models[[2]], i, 1)$MaxAffinity)
  input.models[[i]]$NB[1:4] = input.models[[i]]$NB[1:4] - energy.adjustment
}

# evaluate all sequences and create plots
genome.profile = function(sequence) {
  info = all.models[[1]]
  for (i in 1:nrow(info)) {
    p = plot.score.genome(genomicSequence = DNAString(sequence), fits = all.models[[2]], index = i, mode = 1)
    print(p + ggtitle(paste0(info$Protein[i], " ", info$Dataset[i])))
  }
}

# Max distance
cost.fun.2 = function(scored.list, base.scores) {
  top.aff.prot1 = exp(max(scored.list[[1]])-max(base.scores[[1]]))
  
  # Important to keep the affinity of the first protein the same
  cost = 0
  if (top.aff.prot1>=1) {
    cost = log10(top.aff.prot1)
  } else {
    cost = log10(top.aff.prot1)*10
  }
  min.score.ratio = Inf
  for (i in 2:length(scored.list)) {
    min.score.ratio = min(min.score.ratio, sum(exp(base.scores[[i]]))/sum(exp(scored.list[[i]])))
  }
  if (min.score.ratio>=1) {
    cost = cost + log10(min.score.ratio)
  } else {
    cost = cost + log10(min.score.ratio)*10
  }
  return(max(cost, 0)*10)
}

# Max distance + bonus for significant improvement
cost.fun.3 = function(scored.list, base.scores) {
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
    cost = -log2(-cost)
  }
  return(cost)
}

# simple cost function: Maintain/improve the affinity of the first protein, drop the affinity of the rest
cost.fun = function(scored.list, base.scores) {
  cost = log(sum(exp(scored.list[[1]]))/sum(exp(base.scores[[1]])))
  for (i in 2:length(scored.list)) {
    cost = cost + log(sum(exp(base.scores[[i]]))/sum(exp(scored.list[[i]])))
  }
  return(cost)
}

# Optimize the fkh250 sequence
fkh250 = "CCTCGTCCCACAGCTGGCGATTAATCTTGACATTGAG"
# Loop over all proteins
output.profiles = list()
# for (i in 1:length(input.models)) {
#   # Select a single model
#   idx = 1:length(input.models)
#   idx = c(i, idx[-i])
#   curr.list = input.models[idx]
#   cat(paste0("Optimizing ", names(curr.list)[1], ":\n"))
#   output.profiles[[i]] = optimal.sequence(model.list = curr.list, cost.function = cost.fun, seed.seq = fkh250, verbose = TRUE)
# }

i=16
  # Select a single model
  idx = 1:length(input.models)
  idx = c(i, idx[-i])
  curr.list = input.models[idx]
  cat(paste0("Optimizing ", names(curr.list)[1], ":\n"))
  out3 = optimal.sequence(model.list = curr.list, cost.function = cost.fun.3, seed.seq = fkh250, verbose = TRUE)
  out2 = optimal.sequence(model.list = curr.list, cost.function = cost.fun.2, seed.seq = fkh250, verbose = TRUE)
  out1 = optimal.sequence(model.list = curr.list, cost.function = cost.fun, seed.seq = fkh250, verbose = TRUE)