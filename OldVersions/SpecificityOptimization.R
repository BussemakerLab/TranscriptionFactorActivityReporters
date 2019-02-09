# Define functions and conversions
#Base to number
a.to.n.conversion = c(1:4)
names(a.to.n.conversion) = c("A","C","G","T")
n.to.a.conversion = c("A","C","G","T")

maxSeq = function(k, nuc, dinuc) {
  #Initialize loop (position 1)
  max.list = nuc[,1] #(A C G T)
  char.list= c("A", "C", "G", "T")
  temp.list= char.list
  curr.list = matrix(data=0, 4, 4)
  #Loop over all positions
  for (currPos in 2:k) {
    #Loop over all previous bases
    for (prevBase in 1:4) {
      curr.list[,prevBase] = max.list[prevBase] + nuc[, currPos]
      curr.list[,prevBase] = curr.list[,prevBase] + dinuc[((prevBase-1)*4+1):(prevBase*4), (currPos-1)]
    }
    for (currBase in 1:4) {
      max.list[currBase] = max(curr.list[currBase,])
      temp.list[currBase] = paste0(char.list[which.max(curr.list[currBase,])], c("A", "C", "G", "T")[currBase])
    }
    char.list = temp.list
  }
  return(data.frame(MaxAffinity=max(max.list), BestSeq=char.list[which.max(max.list)], stringsAsFactors = FALSE))
}

# Score returns energy values at each windowed offset  
score = function(seq, nWindows, k, nuc, dinuc) {
  # Create score vector for forward strand and loop over all positions
  score.output = matrix(data=0, nrow=2, ncol=nWindows)
  for (i in 1:nWindows) {
    # extract current windowed sequence, forward and reverse strand
    fwdSeq = seq[i:(i+k-1)]
    revSeq = abs(rev(fwdSeq)-5)
    # loop over all positions in binding site
    fwd = nuc[fwdSeq[1],1]
    rev = nuc[revSeq[1],1]
    for (j in 2:k) {
      fwd = fwd + nuc[fwdSeq[j],j] + dinuc[((fwdSeq[j-1]-1)*4+fwdSeq[j]),j-1]
      rev = rev + nuc[revSeq[j],j] + dinuc[((revSeq[j-1]-1)*4+revSeq[j]),j-1]
    }
    score.output[1,i] = fwd
    score.output[2,i] = rev
  }
  return(score.output)
}

# Optimizers: random-axis optimization
random.axis = function(windows, ks, nuc.vals, dinuc.vals, energy.offsets, start.seq, 
                       cost.function, mutate.dinucs, verbose, mutation.penalty) {
  nModels = length(windows)
  scores = vector(mode="list", length = nModels)
  if (mutate.dinucs) {
    costs = matrix(data = 0,nrow=4, ncol=4)
  } else {
    costs = rep(0, 4)
  }
  
  # Convert sequence to base 4 encoding
  seq.len = length(start.seq)
  curr.opt.seq = rep(0, seq.len)
  for (i in 1:seq.len) {
    curr.opt.seq[i] = a.to.n.conversion[start.seq[i]]
  }
  base.seq = curr.opt.seq
  
  # Get starting value
  for (i in 1:nModels) {
    scores[[i]] = score(curr.opt.seq, windows[i], ks[i], nuc.vals[[i]], dinuc.vals[[i]]) - energy.offsets[i]
  }
  # print(scores)
  # if (TRUE) {
  #   print(start.seq)
  #   print(cost.function(scores))
  #   return();
  # }
  start.cost = cost.function(scores, 0)
  curr.opt.cost = start.cost
  if (verbose) {
    print(paste0("Starting Cost: ", curr.opt.cost))
  }
  
  # Loop over positions randomly until stable solution has been reached
  while (TRUE) {
    last.opt.cost = curr.opt.cost
    if (mutate.dinucs) {
      # Mutate positions randomly with dinucleotide mutations
      mutate.pos = sample(1:(seq.len-1), seq.len-1, replace=F)
      for (curr.pair in mutate.pos) {
        test.seq = curr.opt.seq
        # Test all bases at position
        for (curr.base.A in 1:4) {
          for (curr.base.B in 1:4) {
            test.seq[c(curr.pair, curr.pair+1)] = c(curr.base.A, curr.base.B)
            for (i in 1:nModels) {
              scores[[i]] = score(test.seq, windows[i], ks[i], nuc.vals[[i]], dinuc.vals[[i]]) - energy.offsets[i]
            }
            costs[curr.base.A,curr.base.B] = cost.function(scores, mutation.penalty*sum(test.seq!=base.seq))
          }
        }
        # Update optimal sequence only if the cost has improved
        if (curr.opt.cost < max(costs)) {
          curr.opt.cost = max(costs)
          curr.opt.seq[c(curr.pair, curr.pair+1)] = which(costs == max(costs), arr.ind = TRUE)[1,]
        }
      }
    } else {
      # Mutate positions randomly with nucleotide mutations
      mutate.pos = sample(1:seq.len, seq.len, replace=F)
      for (curr.pos in mutate.pos) {
        test.seq = curr.opt.seq
        # Test all bases at position
        for (curr.base in 1:4) {
          test.seq[curr.pos] = curr.base
          for (i in 1:nModels) {
            scores[[i]] = score(test.seq, windows[i], ks[i], nuc.vals[[i]], dinuc.vals[[i]]) - energy.offsets[i]
          }
          costs[curr.base] = cost.function(scores, mutation.penalty*sum(test.seq!=base.seq))
        }
        # Update optimal sequence only if the cost has improved
        if (curr.opt.cost < max(costs)) {
          curr.opt.cost = max(costs)
          curr.opt.seq[curr.pos] = which.max(costs)
        }
      }
    }
    if (curr.opt.cost == last.opt.cost) {
      break
    }
    if (verbose) {
      print(paste0("Current Cost: ", curr.opt.cost))
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

# Optimizers: simualted annealing
sa = function(windows, ks, nuc.vals, dinuc.vals, energy.offsets, start.seq, cost.function, verbose, mutation.penalty) {
  init.steps = 1000
  init.temp = .9
  temp.decrease = .8
  step.increase = 1.4
  Tcrit = .25
  
  nModels = length(windows)
  scores = vector(mode="list", length = nModels)
  costs = rep(0,4)
  
  # Convert sequence to base 4 encoding
  seq.len = length(start.seq)
  curr.seq = rep(0, seq.len)
  for (i in 1:seq.len) {
    curr.seq[i] = a.to.n.conversion[start.seq[i]]
  }
  base.seq = curr.seq
  
  # Get starting value
  for (i in 1:nModels) {
    scores[[i]] = score(curr.seq, windows[i], ks[i], nuc.vals[[i]], dinuc.vals[[i]]) - energy.offsets[i]
  }
  start.cost = -cost.function(scores, 0)
  curr.cost = start.cost
  if (verbose) {
    cat(paste0("Starting Cost: ", -curr.cost,"\n"))
  }
  
  # Loop over positions randomly until stable solution has been reached
  steps = init.steps
  temp = init.temp
  prev.seq = curr.seq
  prev.cost = curr.cost
  while (temp > Tcrit) {
    # Single epoch
    for (i in 1:steps) {
      # Define next mutation; first decide if two bases or just one is to be mutated
      while (TRUE) {
        # Randomly alternate between mutating one base or two bases
        if (runif(1)>.5) {
          # Two bases
          mutate.pos = sample(1:seq.len, 2, replace=F)
          mutate.base = sample(1:4, 2, replace=T)
          curr.seq[mutate.pos] = mutate.base
        } else {
          # Only one base
          curr.seq[sample(1:seq.len, 1)] = sample(1:4, 1)
        }
        # Accept new move if a different sequence has been generated
        if (sum(curr.seq!=prev.seq)>0) {
          break
        }
      }
      # Compute cost and accept/reject move
      for (i in 1:nModels) {
        scores[[i]] = score(curr.seq, windows[i], ks[i], nuc.vals[[i]], dinuc.vals[[i]]) - energy.offsets[i]
      }
      curr.cost = -cost.function(scores, mutation.penalty*sum(curr.seq!=base.seq))
      alpha = min(1, exp((prev.cost-curr.cost)/temp))
      if ( (curr.cost<prev.cost) || (runif(1)<alpha) ) {
        # Accept move
        prev.seq = curr.seq
        prev.cost = curr.cost
      } else {
        # Reject move
        curr.seq = prev.seq
        curr.cost = prev.cost
      }
    }
    steps = ceiling(steps*step.increase)
    temp = temp*temp.decrease
    if (verbose) {
      cat(paste0("Epoch complete. New T: ",temp,"; New Steps: ",steps,"; Current Cost: ", -curr.cost, "; Curr Seq: ", paste(curr.seq, collapse=""), "\n"))
    }
  }
  
  # Convert optimal sequence to string representation
  opt.seq = start.seq
  for (i in 1:seq.len) {
    opt.seq[i] = n.to.a.conversion[curr.seq[i]]
  }
  opt.seq = paste0(opt.seq, collapse="")
  return(c(-start.cost, opt.seq, -curr.cost))
}

# Sequence optimizer
optimal.sequence = function(..., model.list = NA, cost.function, seq.len = NA, seed.seq = NULL, use.pattern.search = TRUE, 
                            return.max.only = TRUE, verbose = FALSE, penalize.mutations = FALSE) {
  warm.start.multiple = 4
  mutate.dinucs = TRUE

  if (is.na(model.list)) {
    input.models = list(...)
  } else {
    input.models = model.list
  }
  nModels = length(input.models)
  # First, transform fits into vectors, find ks, and the theoretical max relative affinity
  nuc.vals = vector(mode = "list", length = nModels)
  dinuc.vals = vector(mode = "list", length = nModels)
  ks = rep(0, nModels)
  energy.offsets = ks
  windows = ks
  for (i in 1:nModels) {
    nuc = input.models[[i]]$NB
    dinuc = input.models[[i]]$DB
    ks[i] = length(nuc)/4
    dim(nuc) = c(4, ks[i])
    rownames(nuc) = c("A","C","G","T")
    dim(dinuc) = c(16, ks[i]-1)
    rownames(dinuc) = c("AA", "AC", "AG", "AT", 
                        "CA", "CC", "CG", "CT", "GA", "GC", "GG", 
                        "GT", "TA", "TC", "TG", "TT")
    energy.offsets[i] = maxSeq(ks[i], nuc, dinuc)$MaxAffinity
    nuc.vals[[i]] = nuc
    dinuc.vals[[i]] = dinuc
  }
  # Set mutation penalty
  if (!penalize.mutations) {
    mutation.penalty = 0
  } else {
    mutation.penalty = 1
  }
  
  # Preparing the seed sequence depends on the optimizer used

  # Next, optimize the sequence
  # Optimizers: Wang-Landau DOS-style? Typical pattern search (with warm starts)
  
  # First, appropriately handle the length of the string to be optimized and the seed. 
  # This is important due to the use of warm starts. We do so by finding the length of
  # the random 'fields'. Start by handling the case of no seed sequence.
  if (is.null(seed.seq) || is.na(seed.seq)) {
    is.seeded = FALSE
    # Case seq.len to NA if NULL
    if (is.null(seq.len)) {
      seq.len = NA
    }
    # The maximum sequence length is then constrained by the length of the models
    seq.len = max(ks, seq.len, na.rm = TRUE)
    # Random field length is irrelevant
    l.len = seq.len
    r.len = 0
  } else {
    is.seeded = TRUE
    # Max seq length must now consider the length of the seed sequence
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
  }
  # Given the sequence length to optimize over, precompute the windows for each model
  for (i in 1:nModels) {
    windows[i] = seq.len-ks[i]+1;
  }
  if (verbose) {
    print(paste0("Sequence length to optimize: ",seq.len,"; Left Random Field: ", l.len,", Right Random Field: ", r.len))
  }
  
  # Use pattern search (greedy method)
  if (use.pattern.search) {
    # Compute the total number of warm starts to be used (based on the number of variable positions) and loop until completion
    run.outputs = NULL
    for (i in 1:(1 + (l.len+r.len)*warm.start.multiple)) {
      # Construct starting sequence
      start.seq = c(sample(c("A","C","G","T"), l.len, replace = TRUE), seed.seq, sample(c("A","C","G","T"), r.len, replace = TRUE))
      # Run optimization
      output = random.axis(windows, ks, nuc.vals, dinuc.vals, energy.offsets, start.seq, 
                           cost.function, mutate.dinucs, verbose, mutation.penalty)
      # Collect output
      run.outputs = rbind(run.outputs, c(paste(start.seq, collapse=""), output))
    }
    run.outputs = data.frame(run.outputs, stringsAsFactors = FALSE)
    run.outputs[,2] = as.numeric(run.outputs[,2])
    run.outputs[,4] = as.numeric(run.outputs[,4])
    names(run.outputs) = c("StartSequence", "StartCost", "EndSequence", "EndCost")
  } else {
    # Use simulated annealing
    start.seq = c(rep("A", l.len), seed.seq, rep("A", r.len))
    run.outputs = c(paste(start.seq, collapse=""), sa(windows, ks, nuc.vals, dinuc.vals, energy.offsets, start.seq, 
                                                      cost.function, verbose, mutation.penalty))
  }
  if (return.max.only) {
    return(run.outputs[which(run.outputs$EndCost==max(run.outputs$EndCost)),])
  } else {
    return(run.outputs)
  }
}

# Convenient barrier function. Height determines strength of penality, sensitivity is how 'square' the 
# logistic response is, and width is how 'wide' the non-asymptotic regions are
barrier.func = function(height, sensitivity, width, cost) {
  # log(1/99) = -4.59512
  return(height/(1+exp(-4.59512*sensitivity/width*cost)) )
}