wang.landau = function(windows, ks, nuc.vals, dinuc.vals, start.seq, cost.function, verbose, mutation.penalty) {
  # Size of energy divisions
  energy.div.size = 1
  f = 1
  minCount = 10
  fCrit = .8
  tol = .5
  
  nModels = length(windows)
  scores = vector(mode="list", length=nModels)

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
  for (i in 1:nModels) {
    scores[[i]] = score(curr.opt.seq, windows[i], ks[i], nuc.vals[[i]], dinuc.vals[[i]])
  }
  base.scores = scores
  start.cost = cost.function(scores, base.scores, 0)
  En = start.cost
  curr.opt.cost = En
  if (verbose) {
    cat(paste0("Starting Cost: ", start.cost, "\n"))
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
   xn[mutate.pos:(mutate.pos+1)] = sample(1:4, 2, replace=T)
   # nMuts = sample(max.mutate.bases, 1)
   # mutate.pos = sample(seq.pos, nMuts)
   # xn[mutate.pos] = sample(1:4, nMuts, replace=T)
    # Compute new cost
    for (i in 1:nModels) {
      scores[[i]] = score(xn, windows[i], ks[i], nuc.vals[[i]], dinuc.vals[[i]])
    }
    En = cost.function(scores, base.scores, mutation.distance = mutation.penalty*sum(test.seq!=base.seq))
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
#        cat(paste0("; Total Iterations: ", iters, "; mean: ", mean(H[H>0]), "; sd: ", sd(H[H>0]), "; min: ", minH, "\n"))
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
        if (f<=tol) {
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