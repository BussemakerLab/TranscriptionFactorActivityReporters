library(NRLBtools)
library(RapidScore)
library(MultinomialMethods)
# Helper functions
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

# Load models. First find csv
load("~/Documents/GitWorkspaces/NRLB/R-package/NRLBtools/R/sysdata.rda")

# ATF4 model is index 56, CEBP is 57, hetero is 60
# Create rapid score objects
nsb.floor = 2E-4
nSamples = 1E6

ATF4 = new(RapidScore, log(nsb.floor))
CEBP = new(RapidScore, log(nsb.floor))
Het = new(RapidScore, log(nsb.floor))

#ATF4
offset = log(max.seq(NRLBModels, 56)$MaxAffinity[1])
NB = NRLBModels[[2]][[56]]$NB
DB = NRLBModels[[2]][[56]]$DB
NB[1:4] = NB[1:4] - offset
ATF4$add(NB, DB)

#CEBP
offset = log(max.seq(NRLBModels, 57)$MaxAffinity[1])
NB = NRLBModels[[2]][[57]][[1]]$NB
DB = NRLBModels[[2]][[57]][[1]]$DB
NB[1:4] = NB[1:4] - offset
CEBP$add(NB, DB)

#Het
offset = log(max.seq(NRLBModels, 60)$MaxAffinity[1])
NB = NRLBModels[[2]][[60]][[1]]$NB
DB = NRLBModels[[2]][[60]][[1]]$DB
NB[1:4] = NB[1:4] - offset
Het$add(NB, DB)

# Store base sequences
lFlank = convertSequence("GAGTTCTACAGTCCGACGATCCCTGGCGAA")
rFlank = convertSequence("TTCGCCAGGCGTATGT")
ATF4.base.seq = convertSequence("TGACGTCA")
CEBP.base.seq = convertSequence("TTGCGCAA")
Het.base.seq = convertSequence("TGACGCAA")

# Create energy bins and store sequences [Bins are on log-sum-affinity]
model = ATF4
base.seq = ATF4.base.seq

min.val = log(2*(length(c(lFlank, base.seq, rFlank))-length(NB)/4+1)*nsb.floor) # Assumes fixed window size 
seq.len = length(base.seq)
# Create storage matrices
seqs = matrix(data = "NA", nrow = abs(floor(min.val*10)), ncol=(seq.len+1))
affs = matrix(data = NA, nrow = abs(floor(min.val*10)), ncol=(seq.len+1))
# Loop for fixed number of samples
xCurr = base.seq
xPrev = base.seq
base.val = log(sum(exp(model$scoreBulk(c(lFlank, xCurr, rFlank))[[1]])))
prev.val = base.val
for (i in 1:nSamples) {
  # Create new random sequence with a single base mutation
  xCurr[sample(1:seq.len, size = 1)] = sample(0:3, 1)
  tot.ddG = log(sum(exp(model$scoreBulk(c(lFlank, xCurr, rFlank))[[1]])))
  mut.dist = min(sum(xCurr!=base.seq), sum(rev(-1*(xCurr-3))!=base.seq))+1
  aff.dist = abs(floor(tot.ddG*10))
  if (is.na(affs[aff.dist, mut.dist])) {
    # Store newly discovered sequence in bins 
    affs[aff.dist, mut.dist] = tot.ddG
    seqs[aff.dist, mut.dist] = convertSequence(xCurr)
  }
  # MH acceptance step
  if (exp(tot.ddG-prev.val)>runif(1)) {
    # Accept
    xPrev = xCurr
    prev.val = tot.ddG
  } else {
    xCurr = xPrev
  }
}

# Create combined output
ATF4.output = seqs
for (i in 1:nrow(seqs)) {
  for (j in 1:ncol(seqs)) {
    if (is.na(affs[i,j])) {
      next
    }
    ATF4.output[i,j] = paste0(seqs[i,j], " (", round(exp(affs[i,j]-base.val), 4), ")")
  }
}


# Now do CEBP
model = CEBP
base.seq = CEBP.base.seq

min.val = log(2*(length(c(lFlank, base.seq, rFlank))-length(NB)/4+1)*nsb.floor) # Assumes fixed window size 
seq.len = length(base.seq)
# Create storage matrices
seqs = matrix(data = "NA", nrow = abs(floor(min.val*10)), ncol=(seq.len+1))
affs = matrix(data = NA, nrow = abs(floor(min.val*10)), ncol=(seq.len+1))
# Loop for fixed number of samples
xCurr = base.seq
xPrev = base.seq
base.val = log(sum(exp(model$scoreBulk(c(lFlank, xCurr, rFlank))[[1]])))
prev.val = base.val
for (i in 1:nSamples) {
  # Create new random sequence with a single base mutation
  xCurr[sample(1:seq.len, size = 1)] = sample(0:3, 1)
  tot.ddG = log(sum(exp(model$scoreBulk(c(lFlank, xCurr, rFlank))[[1]])))
  mut.dist = min(sum(xCurr!=base.seq), sum(rev(-1*(xCurr-3))!=base.seq))+1
  aff.dist = abs(floor(tot.ddG*10))
  if (is.na(affs[aff.dist, mut.dist])) {
    # Store newly discovered sequence in bins 
    affs[aff.dist, mut.dist] = tot.ddG
    seqs[aff.dist, mut.dist] = convertSequence(xCurr)
  }
  # MH acceptance step
  if (exp(tot.ddG-prev.val)>runif(1)) {
    # Accept
    xPrev = xCurr
    prev.val = tot.ddG
  } else {
    xCurr = xPrev
  }
}

# Create combined output
CEBP.output = seqs
for (i in 1:nrow(seqs)) {
  for (j in 1:ncol(seqs)) {
    if (is.na(affs[i,j])) {
      next
    }
    CEBP.output[i,j] = paste0(seqs[i,j], " (", round(exp(affs[i,j]-base.val), 4), ")")
  }
}


# Now do CEBP
model = Het
base.seq = Het.base.seq

min.val = log(2*(length(c(lFlank, base.seq, rFlank))-length(NB)/4+1)*nsb.floor) # Assumes fixed window size 
seq.len = length(base.seq)
# Create storage matrices
seqs = matrix(data = "NA", nrow = abs(floor(min.val*10)), ncol=(seq.len+1))
affs = matrix(data = NA, nrow = abs(floor(min.val*10)), ncol=(seq.len+1))
# Loop for fixed number of samples
xCurr = base.seq
xPrev = base.seq
base.val = log(sum(exp(model$scoreBulk(c(lFlank, xCurr, rFlank))[[1]])))
prev.val = base.val
for (i in 1:nSamples) {
  # Create new random sequence with a single base mutation
  xCurr[sample(1:seq.len, size = 1)] = sample(0:3, 1)
  tot.ddG = log(sum(exp(model$scoreBulk(c(lFlank, xCurr, rFlank))[[1]])))
  mut.dist = min(sum(xCurr!=base.seq), sum(rev(-1*(xCurr-3))!=base.seq))+1
  aff.dist = abs(floor(tot.ddG*10))
  if (is.na(affs[aff.dist, mut.dist])) {
    # Store newly discovered sequence in bins 
    affs[aff.dist, mut.dist] = tot.ddG
    seqs[aff.dist, mut.dist] = convertSequence(xCurr)
  }
  # MH acceptance step
  if (exp(tot.ddG-prev.val)>runif(1)) {
    # Accept
    xPrev = xCurr
    prev.val = tot.ddG
  } else {
    xCurr = xPrev
  }
}

# Create combined output
Het.output = seqs
for (i in 1:nrow(seqs)) {
  for (j in 1:ncol(seqs)) {
    if (is.na(affs[i,j])) {
      next
    }
    Het.output[i,j] = paste0(seqs[i,j], " (", round(exp(affs[i,j]-base.val), 4), ")")
  }
}

write.table(x = ATF4.output, file = "~/Desktop/ATF4.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = CEBP.output, file = "~/Desktop/CEBP.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = Het.output, file = "~/Desktop/Het.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
