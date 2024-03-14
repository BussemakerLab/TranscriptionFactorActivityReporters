setwd("~/Documents/Research/SequenceEngineering/NKI/")
load("v1Run.rda")
source("~/Documents/Research/SequenceOptimization/base/WLDOS_SequenceOptimizer.R")
source("~/Documents/Research/SequenceOptimization/base/PlotLibrary.R")
load("/Users/chaitanya/Documents/GitWorkspaces/NRLB/R-package/NRLBtools/R/sysdata.rda")
library(NRLBtools)
library(jsonlite)
library(MultinomialMethods)

# Load p53
input.models = list()
idx = which(NRLBModelInfo$Protein=="p53" & NRLBModelInfo$Info=="WildType")
input.models[[1]] = NRLBModels[[2]][[idx]][[1]]
energy.adjustment = log(max.seq(NRLBModels, idx, 1)$MaxAffinity)
input.models[[1]]$NB[1:4] = input.models[[1]]$NB[1:4] - energy.adjustment
names(input.models) = "p53"
opt.seq = max.seq(NRLBModels, idx, 1)$BestSeq

seqs = p53.cost.fun.out$SequenceExamples
seqs[seqs=="NA"] = NA
vals = score.seq.helper(sequences = seqs, base.seq = opt.seq, input.models)
vals = log(vals$Sum)
df = melt(vals, value.name="Score", varnames = c("Mutations", "Protein"))
df$Mutations = rep(1:nrow(seqs)-1, ncol(seqs))
df$ddG = rep(rev(-1*(1:ncol(seqs)-1)), each=nrow(seqs))

# Split by ddG
df = split(df, df$ddG)
output = NULL
for (i in 1:length(df)) {
  row.names(df[[i]]) = NULL
  curr.val = df[[i]][which.min(round(abs(df[[i]]$Score-df[[i]]$ddG)*100)+df[[i]]$Mutations*2),]
  output = rbind(output, cbind(curr.val, Sequence=seqs[as.numeric(row.names(curr.val)),i]))
}


# Load GR
json.con = file(description = "fit.Zhang2017.GR.7297.json", open="r")
GR.json = jsonlite::stream_in(json.con)
close(json.con)
# Extract betas
GR.bm = GR.json$coefficients$bindingModes[[1]]
mono = GR.bm$mononucleotide[[2]]
di = as.numeric(GR.bm$dinucleotide[[2]])
# load into previous input.models
input.models[[1]]$NB = mono
input.models[[1]]$DB = di
names(input.models) = "GR"
NRLBModels[[2]][[idx]][[1]]$NB = mono
NRLBModels[[2]][[idx]][[1]]$DB = di
energy.adjustment = log(max.seq(NRLBModels, idx, 1)$MaxAffinity)
input.models[[1]]$NB[1:4] = input.models[[1]]$NB[1:4] - energy.adjustment
opt.seq = max.seq(NRLBModels, idx, 1)$BestSeq

seqs = GR.cost.fun.out$SequenceExamples
seqs[seqs=="NA"] = NA
seqs = seqs[,2:(ncol(seqs)-1)]
vals = score.seq.helper(sequences = seqs, base.seq = opt.seq, input.models)
vals = log(vals$Sum)
df = melt(vals, value.name="Score", varnames = c("Mutations", "Protein"))
df$Mutations = rep(1:nrow(seqs)-1, ncol(seqs))
df$ddG = rep(rev(-.5*(1:ncol(seqs)-1)), each=nrow(seqs))

# Split by ddG
df = split(df, df$ddG)
for (i in 1:length(df)) {
  row.names(df[[i]]) = NULL
  curr.val = df[[i]][which.min(round(abs(df[[i]]$Score-df[[i]]$ddG)*100)+df[[i]]$Mutations*2),]
  output = rbind(output, cbind(curr.val, Sequence=seqs[as.numeric(row.names(curr.val)),i]))
}

output = output[order(output$Protein, -output$ddG),]
row.names(output) = NULL
