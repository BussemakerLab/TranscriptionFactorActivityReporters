# Must install RapidScore C++ Library First!

library(RapidScore)
library(abind)
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

# biased 2D wang-landau (looks only for positive mutations) - REMOVE
wang.landau = function(score.object, start.seq, cost.function, isBiased=TRUE, verbose=TRUE) {
  verbose = FALSE
  isBiased = TRUE
  # Size of energy divisions and other parameters - play with these?
  energy.div.size = 1
  f = 1
  minCount = 100
  minIters = 1E4
  fCrit = .25
  tol = .5
  
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
  # State variables
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
  prevHIdx = NA
  prevOptCost = NA
  
  # Initialize Wang-Landau Scheme
  IE0= round((E0-Emin)/energy.div.size)+1
  IM0= 1  # Initial bin for edit distance is 0 
  S[IM0, IE0] = S[IM0, IE0] + f
  H[IM0, IE0] = H[IM0, IE0] + 1
  HIdx[IM0, IE0] = TRUE
  last.seq[IM0, IE0, ] = x0
  iters = 1
  
  # Loop continuously until convergence
  while (TRUE) {
    # Propose a new state with a nucleotide mutation. On special steps, move to a 
    # bin with fewer visits (hopefully). This biased move needs to be compensated
    if (runif(1) <= .005) {   # TODO: does this rate need to scale with sequence length?
      # Compute bin sampling bias distributon
      cts = c(H)
      cts[cts>0] = sum(cts)/cts[cts>0]
      cts = cts/sum(cts)
      idx = arrayInd(sample(length(H), 1, prob = cts), dim(H))
      xn = last.seq[idx[1],idx[2],]
      # This transition ratio is to account for unbalanced transition probabilities 
      trans.ratio = cts[(IE0-1)*mDivs+IM0]/cts[(idx[2]-1)*mDivs+idx[1]] 
    } else {
      # 'Standard' move where the sequence is mutated at a single position
      xn = x0
      mutate.pos = sample(seq.pos, 1)
      xn[mutate.pos] = sample(0:3, 1)
      trans.ratio = 1
    }
    
    # Compute new cost
    scores = score.object$scoreBulk(xn)
    En = cost.function(scores, base.scores)
    # Biased sampling: Reject if sequence energy is less than the least value TODO: REMOVE!!!
    if (isBiased && En < Emin) {
      next
    }
    # Compute mutation distance; Reverse complement symmetry must be accounted for when computing
    Mn = min(sum(xn!=base.seq), sum(rev(-1*(xn-3))!=base.seq)) 
    # Compute new cost and mutation indexes
    IEn = round((En-Emin)/energy.div.size)+1
    IMn = Mn+1
    
    # Update ranges and restart wang-landau process if a sequence is produced that is out of range
    m.range.update = FALSE
    e.range.update = FALSE
    # Check mutation ranges
    if (IMn > mDivs) {
      m.range.update = TRUE
      if (verbose) {
        cat(paste0("Old Max Mutation: ", MutMax, "; New Max Mutation: ", IMn-1, "; "))
      }
      # Update and add rows to opt.seq
      for (k in (MutMax+1):(IMn-1)) {
        opt.seqs = rbind(opt.seqs, rep(c(k, -Inf, list(NA))))
      }
      # Recompute Range and divs
      MutMax = IMn-1
      mDivs = IMn
    }
    # Check cost ranges
    if (IEn<=0 || IEn > eDivs) {
      e.range.update = TRUE
      if (verbose) {
        cat(paste0("Old Cost Range: ", Emin, ":", Emax,"; "))
      }
      # Recompute range and divs
      Emin.new = floor(min(Emin, En)/energy.div.size)*energy.div.size
      Emax.new = ceiling(max(Emax, En)/energy.div.size)*energy.div.size
      if (verbose) {
        cat(paste0("New Cost Range: ", Emin.new, ":", Emax.new))
      }
    }
    # Update range and restart wang-landau process if necessary
    if (m.range.update || e.range.update) {
      # Do not need to restart
      if (f==1) {
        oldHIdx = sum(HIdx)
        minH = min(H[H>0])
        # Begin by padding the H and S matrices with mutation range increases
        if (m.range.update) {
          last.seq = abind(last.seq, array(dim=c(mDivs-nrow(H), ncol(H), seq.len)), along = 1)
          H = rbind(H, matrix(data = 0, ncol = ncol(H), nrow = mDivs-nrow(H)))
          S = rbind(S, matrix(data = 0, ncol = ncol(S), nrow = mDivs-nrow(S)))
          HIdx = rbind(HIdx, matrix(data = FALSE, ncol = ncol(HIdx), nrow = mDivs-nrow(HIdx)))
        }
        # Pad H and S with energy range increases
        if (e.range.update) {
          H = cbind(matrix(data = 0, ncol = Emin-Emin.new, nrow = nrow(H)), H, matrix(data = 0, ncol = Emax.new-Emax, nrow = nrow(H)))
          S = cbind(matrix(data = 0, ncol = Emin-Emin.new, nrow = nrow(S)), S, matrix(data = 0, ncol = Emax.new-Emax, nrow = nrow(S)))
          HIdx = cbind(matrix(data = FALSE, ncol = Emin-Emin.new, nrow = nrow(HIdx)), HIdx, matrix(data = FALSE, ncol = Emax.new-Emax, nrow = nrow(HIdx)))
          last.seq = abind(array(dim=c(nrow(H), Emin-Emin.new, seq.len)), last.seq, array(dim=c(nrow(H), Emax.new-Emax, seq.len)), along=2)
          eDivs = ncol(H)
          Emin = Emin.new
          Emax = Emax.new
          IEn = round((En-Emin)/energy.div.size)+1
        }
        if (verbose) {
          cat(paste0("; Total Iterations: ", iters, "; mean: ", mean(H[H>0]), "; sd: ", sd(H[H>0]), "; min: ", minH, "; oldHIdx: ", oldHIdx, "; newHIdx: ", sum(HIdx), "\n"))
        }
      } else {
        # Need a hard restart, BUT PRESERVE HIdx
        mDivs = MutMax+1
        HIdx = HIdx | prevHIdx
        if (m.range.update) {
          HIdx = rbind(HIdx, matrix(data = FALSE, ncol = ncol(HIdx), nrow = mDivs-nrow(HIdx)))
        }
        if (e.range.update) {
          HIdx = cbind(matrix(data = FALSE, ncol = Emin-Emin.new, nrow = nrow(HIdx)), HIdx, matrix(data = FALSE, ncol = Emax.new-Emax, nrow = nrow(HIdx)))
          Emin = Emin.new
          Emax = Emax.new
          eDivs =((Emax-Emin)/energy.div.size+1)
        }
        # This order of operations preserves active bins, both in current and previous steps
        H = matrix(data = 0, ncol = eDivs, nrow = mDivs)
        S = H
        x0= xn
        E0= En
        IE0= round((E0-Emin)/energy.div.size)+1
        IM0= IMn
        iters = 0
        f = 1
        S[IM0, IE0] = S[IM0, IE0] + f
        H[IM0, IE0] = H[IM0, IE0] + 1
        HIdx[IM0, IE0] = TRUE
        last.seq = array(dim = c(mDivs, eDivs, seq.len))
        last.seq[IM0, IE0, ] = x0
        cat(paste0("prevHIdx: ", sum(prevHIdx), "; HIdx: ", sum(HIdx), ": "))
        prevHIdx = NA
        cat("Energy Ranges Updated - Restarted sampling!\n")
        next
      }
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
    }
    # Update
    S[IM0, IE0] = S[IM0, IE0] + f
    H[IM0, IE0] = H[IM0, IE0] + 1
    iters = iters + 1
    
    # Monitoring
    if (iters %% 250000 ==0) {
      cat(paste0("Total Iterations: ", iters, "; mean: ", mean(H[HIdx]), "; sd: ", sd(H[HIdx]), "; min: ", min(H[HIdx]), "; active bins: ", sum(HIdx), "; # Zero Bins: ", sum(H[HIdx]==0), "\n"))
      # If the number of active bins changes after f goes below 1, need to restart sampling
      if (f<1 && sum(!(which(HIdx) %in% which(prevHIdx)))>0) {
        # NOTE: The ranges of H and S DO NOT change here; its that a new 'bin' has opened up 
        # and sampling needs to begin again; begin by resetting state variables
        f = 1
        H = matrix(data = 0, ncol = eDivs, nrow = mDivs)
        S = H
        # HIdx now needs to contain a merged value - This is important! Ensures that all bins are tracked
        cat(paste0("prevHIdx: ", sum(prevHIdx), "; HIdx: ", sum(HIdx), " new Idxes: ", sum(!(which(HIdx) %in% which(prevHIdx))), ": "))
        HIdx = HIdx | prevHIdx
        # Reset the last visited seq array
        last.seq = array(dim = c(mDivs, eDivs, seq.len))
        last.seq[IM0, IE0, ] = x0
        # Reset the previous cycle monitors
        prevHIdx = NA
        prevOptCost = NA
        # Re-init
        S[IM0, IE0] = S[IM0, IE0] + f
        H[IM0, IE0] = H[IM0, IE0] + 1
        HIdx[IM0, IE0] = TRUE
        iters = 1
        cat("New active bins detected - Restarted Sampling!\n")
      }
    }
    
    # Ensure that all energy bins have been visited an even number of times
    if (f<1 && min(H[prevHIdx]) > minCount && min(H[prevHIdx]) >= tol*mean(H[prevHIdx]) && iters>minIters) {
      cat(paste0("Flatness criteria reached for f = ",f," in ",iters," iterations. "))
      cat(paste0("Total Active Bins: ", sum(HIdx), "; Previous Round Total Active Bins: ", sum(prevHIdx), "\n"))
      f = f/2;
      # Ensure that the opt cost is stable before exiting
      if (f<=fCrit && !is.na(prevOptCost) && sum(prevOptCost!=opt.seqs$OptimalCost)==0) {
        # Complete, exit
        break
      }
      # Store this round's optimal cost and active bin indices
      prevOptCost = opt.seqs$OptimalCost
      prevHIdx = HIdx
      # Reset the sampling for the next round
      H = matrix(data = 0, ncol = eDivs, nrow = mDivs)
      last.seq = array(dim = c(mDivs, eDivs, seq.len))
      last.seq[IM0, IE0, ] = x0
      S[IM0, IE0] = S[IM0, IE0] + f
      H[IM0, IE0] = H[IM0, IE0] + 1
      HIdx[IM0, IE0] = TRUE
      iters = 1;
    } else if (f==1 && min(H[HIdx]) > minCount && min(H[HIdx]) >= tol*mean(H[HIdx]) && iters>minIters) {
      cat(paste0("Flatness criteria reached for f = 1 in ",iters," iterations. "))
      cat(paste0("Total Active Bins: ", sum(HIdx), "; Previous Round Total Active Bins: ", sum(prevHIdx), "\n"))
      f = f/2;
      # Store this round's optimal cost and active bin indices
      prevOptCost = opt.seqs$OptimalCost
      prevHIdx = HIdx
      # Reset the sampling for the next round
      H = matrix(data = 0, ncol = eDivs, nrow = mDivs)
      last.seq = array(dim = c(mDivs, eDivs, seq.len))
      last.seq[IM0, IE0, ] = x0
      S[IM0, IE0] = S[IM0, IE0] + f
      H[IM0, IE0] = H[IM0, IE0] + 1
      HIdx[IM0, IE0] = TRUE
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
optimal.sequence = function(..., model.list = NA, cost.function, seq.len = NA, seed.seq, verbose = FALSE) {
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
  
  # Optimize using wang-landau sampling. Construct starting sequence
  start.seq = c(sample(c("A","C","G","T"), l.len, replace = TRUE), seed.seq, sample(c("A","C","G","T"), r.len, replace = TRUE))
  # Run optimization
  output = wang.landau(score.object, start.seq, cost.function, isBiased = TRUE, verbose = verbose)
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
cost.fun = function(scored.list, base.scores) {
  cost = log(sum(exp(scored.list[[1]]))/sum(exp(base.scores[[1]]))) - 
    barrier.func(50, 1, exp(max(scored.list[[1]]))/10, 
                 .85*exp(max(base.scores[[1]]))-exp(max(scored.list[[1]])))
  for (i in 2:length(scored.list)) {
    cost = cost + log(sum(exp(base.scores[[i]]))/sum(exp(scored.list[[i]])))
  }
  return(cost)
}

# Optimize the fkh250 sequence
fkh250 = "CCTCGTCCCACAGCTGGCGATTAATCTTGACATTGAG"
# Loop over all proteins
output.profiles = list()
for (i in 1:length(input.models)) {
  # Select a single model
  idx = 1:length(input.models)
  idx = c(i, idx[-i])
  curr.list = input.models[idx]
  cat(paste0("Optimizing ", names(curr.list)[1], ":\n"))
  out = optimal.sequence(model.list = curr.list, cost.function = cost.fun, seed.seq = fkh250, verbose = TRUE)
  output.profiles = c(output.profiles, out)
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