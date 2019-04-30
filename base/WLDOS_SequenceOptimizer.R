# Must install RapidScore C++ Library First!
library(RapidScore)
library(abind)

########################
### Define Functions ###
########################
# Base to number
convertSequence = function(input) {
  a.to.n.conversion = c(0:3)
  names(a.to.n.conversion) = c("A","C","G","T")
  n.to.a.conversion = c("A","C","G","T")
  
  if (class(input)=="character") {
    # Input is text - convert to numeric. See if it needs to be split up
    if (length(input)==1) {
      in.seq = substring(input, 1:nchar(input), 1:nchar(input))
    } else {
      in.seq = input
    }
    seq.len = length(in.seq)
    out.seq = rep(0, seq.len)
    for (i in 1:seq.len) {
      out.seq[i] = a.to.n.conversion[in.seq[i]]
    }
  } else {
    # Input is numeric
    out.seq = input
    for (j in 1:length(input)) {
      out.seq[j] = n.to.a.conversion[input[j]+1]
    }
    out.seq = paste0(out.seq, collapse="")
  }
  return(out.seq)
}

wang.landau = function(score.object, start.seq, cost.function, verbose=TRUE) {
  # Size of energy divisions and other parameters - play with these?
  energy.div.size = 1
  f = 1
  minCount = 100
  minIters = 1E4
  tol = .8
  globalIterations = 0
  nWalkers = 20
  testIterations = 2.5E3
  acceptanceRate = .002
  maxHits = 5
  
  # Convert sequence to base 4 encoding, and define positions to mutate
  seq.len = length(start.seq)
  seq.pos = 1:seq.len
  x0 = matrix(data = rep(convertSequence(start.seq), nWalkers), nrow=nWalkers, byrow = T)
  base.seq = x0[1,]    # Only needs to be a 'scalar"
  xn = x0[1,]
  
  # Compute starting values for all walkers
  scores = score.object$scoreBulk(base.seq)
  base.scores = scores
  start.cost = cost.function(scores, base.scores)
  En = start.cost
  E0 = En
  if (verbose) {
    cat(paste0("Starting Cost: ", start.cost, "\n"))
  }
  
  # Initialize Wang-Landau scheme; start by defining ranges for all wakers
  Emin = floor(min(0, En)/energy.div.size)*energy.div.size
  Emax = ceiling(max(0, En)/energy.div.size)*energy.div.size
  MutMin = 0
  MutMax = 2
  # Define bin size
  eDivs= ((Emax-Emin)/energy.div.size+1)
  mDivs = MutMax+1
  # State variables; HIdx is a binary matrix storing the indicies of visited bins for all walkers
  H = array(dim = c(mDivs, eDivs, nWalkers), data=0)
  S = array(dim = c(mDivs, eDivs, nWalkers), data=0)
  HIdx = matrix(data = FALSE, ncol = eDivs, nrow = mDivs)   # This does not need to 'split up' amongst walkers
  last.seq = array(dim = c(mDivs, eDivs, seq.len))          # Last visited seq array
  
  # List of optimal sequences by mutational distance
  opt.seqs = data.frame(Mutations=MutMin:MutMax, OptimalCost=-Inf)
  opt.seqs$OptimalSeq=list(NA)
  opt.seqs$OptimalCost[1] = E0
  opt.seqs$OptimalSeq[[1]]= x0[1,]
  
  # Monitor previous cycles (is this needed??)
  prevOptCost = NA
  prevHIdx = HIdx
  prevCycleActiveBins = 0
  prevCycleEmptyBins = rep(0, nWalkers)
  prevCycleMinBins   = rep(0, nWalkers)
  hit.counter = 1
  
  # Initialize Wang-Landau Scheme
  IE0= rep(round((E0-Emin)/energy.div.size)+1, nWalkers)
  IM0= rep(1, nWalkers)  # Initial bin for edit distance is 0 
  S[IM0[1], IE0[1],] = S[IM0[1], IE0[1], ] + f
  H[IM0[1], IE0[1],] = H[IM0[1], IE0[1], ] + 1
  HIdx[IM0[1], IE0[1]] = TRUE
  prevHIdx[IM0[1], IE0[1]] = TRUE
  last.seq[IM0[1], IE0[1], ] = x0[1,]
  loops = 0
  iters = 0
  
  # Monitors min, mean, and sd
  walker.monitors = matrix(data=0, nrow=nWalkers, ncol=4)
  colnames(walker.monitors) = c("Empty", "Min", "Mean", "SD")
  
  # Loop continuously until convergence
  while (TRUE) {
    # Loop sequentially through all walkers
    for (currWalker in 1:nWalkers) {
      # Perform fixed number of iterations
      for (v in 1:testIterations) {
        # Propose a new state with a nucleotide mutation. On special steps, move to a 
        # bin with fewer visits (hopefully). This biased move needs to be compensated
        if (hit.counter>=maxHits && runif(1) <= acceptanceRate) {   # TODO: does this rate need to scale with sequence length?
          # Compute bin sampling bias distribution
          Hp = H[ , , currWalker]
          Hp[prevHIdx] = Hp[prevHIdx]-min(Hp[prevHIdx])+1
          cts = c(Hp)
          cts[cts>0] = sum(cts)/cts[cts>0]
          cts = cts/sum(cts)
          idx = arrayInd(sample(length(H[, ,currWalker]), 1, prob = cts), dim(H[ , , currWalker]))
          xn = last.seq[idx[1],idx[2],]
          # xn[sample(seq.pos,1)] = sample(0:3,1) # Only need this in case we are not using parallel tempering
          # This transition ratio is to account for unbalanced transition probabilities
          trans.ratio = cts[(IE0[currWalker]-1)*mDivs+IM0[currWalker]]/cts[(idx[2]-1)*mDivs+idx[1]]
        } else {
          # Propose a new state with a nucleotide mutation
          xn = x0[currWalker, ]
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
        
        # Begin by padding the H and S matrices with mutation range increases
        if (IMn > mDivs) {
          # Update and add rows to opt.seq
          for (k in (MutMax+1):(IMn-1)) {
            opt.seqs = rbind(opt.seqs, rep(c(k, -Inf, list(NA))))
          }
          # Recompute Range and divs
          MutMax = IMn-1
          oldMDivs = mDivs
          mDivs = IMn
          H = abind(H, array(dim=c(mDivs-oldMDivs, eDivs, nWalkers), data = 0), along=1)
          S = abind(S, array(dim=c(mDivs-oldMDivs, eDivs, nWalkers), data = 0), along=1)
          HIdx = rbind(HIdx, matrix(data = FALSE, ncol = eDivs, nrow = mDivs-nrow(HIdx)))
          prevHIdx = rbind(prevHIdx, matrix(data = FALSE, ncol = eDivs, nrow = mDivs-oldMDivs))
          last.seq = abind(last.seq, array(dim=c(mDivs-oldMDivs, eDivs, seq.len), data=0), along = 1)
        }
        # Check cost ranges; pad H and S with energy range increases
        if (IEn<=0 || IEn > eDivs) {
          # Recompute range and divs
          Emin.new = floor(min(Emin, En)/energy.div.size)*energy.div.size
          Emax.new = ceiling(max(Emax, En)/energy.div.size)*energy.div.size
          H = abind(array(dim=c(mDivs, Emin-Emin.new, nWalkers), data=0), H, array(dim=c(mDivs, Emax.new-Emax, nWalkers), data=0), along=2)
          S = abind(array(dim=c(mDivs, Emin-Emin.new, nWalkers), data=0), S, array(dim=c(mDivs, Emax.new-Emax, nWalkers), data=0), along=2)
          HIdx = cbind(matrix(data = FALSE, ncol = Emin-Emin.new, nrow = mDivs), HIdx, matrix(data = FALSE, ncol = Emax.new-Emax, nrow = mDivs))
          prevHIdx = cbind(matrix(data = FALSE, ncol = Emin-Emin.new, nrow = mDivs), 
                           prevHIdx, matrix(data = FALSE, ncol = Emax.new-Emax, nrow = mDivs))
          last.seq = abind(array(dim=c(mDivs, Emin-Emin.new, seq.len), data=0), last.seq, array(dim=c(mDivs, Emax.new-Emax, seq.len), data=0), along=2)
          eDivs = ncol(H[,,1])
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
        a = exp(S[IM0[currWalker], IE0[currWalker], currWalker] - S[IMn, IEn, currWalker])*trans.ratio
        
        # Metropolis-Hastings step
        if (runif(1) <= a) {
          # Accept current proposal
          x0[currWalker,] = xn
          # E0 = En
          IE0[currWalker] = IEn
          IM0[currWalker] = IMn
          HIdx[IMn, IEn] = TRUE
          prevHIdx[IMn, IEn] = TRUE
          last.seq[IMn, IEn, ] = xn
        }
        
        # Update
        S[IM0[currWalker], IE0[currWalker], currWalker] = S[IM0[currWalker], IE0[currWalker], currWalker] + f
        H[IM0[currWalker], IE0[currWalker], currWalker] = H[IM0[currWalker], IE0[currWalker], currWalker] + 1
      }
    }
    
    # Monitor and terminate if all energy bins have been visited an even number of times
    for (currWalker in 1:nWalkers) {
      currH = H[,,currWalker][prevHIdx]
      walker.monitors[currWalker,] = c(sum(currH==0), min(currH[currH>0]), mean(currH[currH>0]), sd(currH[currH>0]))
    }
    iters = iters + testIterations
    globalIterations = globalIterations + testIterations*nWalkers
    
    # Check to see if the number of empty bins changes.
    if (sum(prevHIdx)==prevCycleActiveBins && 
        sum(abs(prevCycleEmptyBins-walker.monitors[,1]))==0 &&
        sum(abs(prevCycleMinBins-walker.monitors[,2]))==0) {
      hit.counter = hit.counter + 1
    } else {
      hit.counter = 1
    }
    prevCycleActiveBins = sum(prevHIdx)
    prevCycleEmptyBins  = walker.monitors[,1]
    prevCycleMinBins    = walker.monitors[,2]
    
    # Print cycle update
    if (globalIterations %% 5E5 == 0) {
      cat(paste0("Total Iterations: ", globalIterations, "; active bins: ", sum(prevHIdx), 
                 "; empty bins: ", sum(HIdx[prevHIdx]==0), "; Emin-Emax: ", Emin, "-", Emax, "; Hit Counter: ", hit.counter,"\n"))
      cat("Number Empty: ")
      print(walker.monitors[,1])
      cat("Min:          ")
      print(walker.monitors[,2])
      # cat("Mean:         ")
      # print(walker.monitors[,3])
      # cat("SD:           ")
      # print(walker.monitors[,4])
    }
    
    if (all(walker.monitors[,"Min"]>minCount) && all(walker.monitors[,"Min"] >= tol*walker.monitors[,"Mean"]) && iters>minIters) {
      # Converged 
      loops = loops + 1
      if (verbose) {
        cat(paste0("Flatness criteria reached in ",iters," individual walker iterations. ", sum(prevHIdx), " active bins; ", loops, " loops completed.\n"))
      }
      print(opt.seqs$OptimalCost)
      # Terminate if the minimum number of loops is met and if the optimal cost vector has remained constant for last loop
      if (loops>1 && sum(prevOptCost!=opt.seqs$OptimalCost)==0) {
        cat(paste0("Finished in ", globalIterations, " total iterations.\n"))
        break
      }
      # Reset the sampling for the next round
      prevOptCost = opt.seqs$OptimalCost
      H = H*0
      S = S*0
      HIdx = H[,,1]>0
      for (currWalker in 1:nWalkers) {
        S[IM0[currWalker], IE0[currWalker], ] = S[IM0[currWalker], IE0[currWalker], ] + f
        H[IM0[currWalker], IE0[currWalker], ] = H[IM0[currWalker], IE0[currWalker], ] + 1
        HIdx[IM0[currWalker], IE0[currWalker]] = TRUE
        prevHIdx[IM0[currWalker], IE0[currWalker]] = TRUE
      }
      iters = 0;
    }
    
    # Swap Chains (Parallel Tempering)
    # Generate exchange pairs
    swap.pairs = matrix(sample(1:nWalkers), ncol=2)
    for (i in 1:nrow(swap.pairs)) {
      # See if pair will be swapped
      x = swap.pairs[i,1]
      y = swap.pairs[i,2]
      a = S[IM0[x], IE0[x], x]*S[IM0[y],IE0[y],y]/(S[IM0[x], IE0[x], y]*S[IM0[y],IE0[y],x])
      if (is.nan(a)) {
        # In case both of the bins have no entropy, always accept (it as NOT been explored in that replicate)
        cat("NaN in replica exchange!\n")
        a = 1
      }
      if (runif(1) <= a) {
        # Swap state
        swapM = IM0[x]
        swapE = IE0[x]
        swapX = x0[x,]
        IM0[x] = IM0[y]
        IE0[x] = IE0[y]
        x0[x,] = x0[y,]
        IM0[y] = swapM
        IE0[y] = swapE
        x0[y,] = swapX
      }
    }
  }
  
  # Convert optimal sequence to string representation
  for (i in 1:nrow(opt.seqs)) {
    curr.seq = opt.seqs$OptimalSeq[[i]]
    if (sum(is.na(curr.seq))>0) {
      next
    }
    opt.seqs$OptimalSeq[[i]] = convertSequence(curr.seq)
  }
  opt.seqs$OptimalSeq = unlist(opt.seqs$OptimalSeq)
  print(opt.seqs)
  
  # Convert last.seq into a matrix
  example.seqs = matrix(data="NA", nrow = nrow(HIdx), ncol = ncol(HIdx))
  for (y in 1:nrow(HIdx)) {
    for (x in 1:ncol(HIdx)) {
      if (HIdx[y,x]) {
        example.seqs[y,x] = convertSequence(last.seq[y,x,])
      }
    }
  }
  return(list(MaximalSequences=opt.seqs, SequenceExamples=example.seqs))
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