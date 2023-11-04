library(NRLBtools)
library(jsonlite)
# Load depending on OS
if (as.character(Sys.info()['sysname'])=="Windows") {
  # Sequence Engineering Library
  source("~/Columbia Work/Research/SequenceOptimization/base/WLDOS_SequenceOptimizer.R")
  # Load p53
  stop("No path to .rda file")
} else {
  # Sequence Engineering Library
  source("~/Documents/Research/SequenceOptimization/base/WLDOS_SequenceOptimizer.R")
  source("~/Documents/Research/SequenceOptimization/base/PlotLibrary.R")
  setwd("~/Documents/Research/SequenceEngineering/NKI/")
  # Load p53
  load("/Users/chaitanya/Documents/GitWorkspaces/NRLB/R-package/NRLBtools/R/sysdata.rda")
}

# Load p53
input.models = list()
idx = which(NRLBModelInfo$Protein=="p53" & NRLBModelInfo$Info=="WildType")
input.models[[1]] = NRLBModels[[2]][[idx]][[1]]
energy.adjustment = log(max.seq(NRLBModels, idx, 1)$MaxAffinity)
input.models[[1]]$NB[1:4] = input.models[[1]]$NB[1:4] - energy.adjustment
names(input.models) = "p53"


# Cost function is input-output mapping
cost.fun = function(scored.list, base.scores) {
  return(scored.list[[1]][1])
}

#########################
### Run Optimizations ###
#########################
# Start with optimal sequence
opt.seq = max.seq(NRLBModels, idx, 1)$BestSeq
p53.cost.fun.out = optimal.sequence(model.list = input.models, cost.function = cost.fun, seed.seq = opt.seq, verbose=TRUE)

plot.frontier(data = p53.cost.fun.out, derivative = F)
ggsave("p53_frontier_v1.pdf", width = 6, height=6)
plot.cost.range(data = p53.cost.fun.out) + xlab(expression(Delta*Delta*G)) + coord_flip()
ggsave("p53_costrange_v1.png", device = "png", width = 6, height=6)

seqs = p53.cost.fun.out$SequenceExamples
seqs[seqs=="NA"] = NA
vals = score.seq.helper(sequences = seqs, base.seq = opt.seq, input.models)
vals = log(vals$Sum)
df = melt(vals, value.name="Score", varnames = c("Mutations", "Protein"))
df$Mutations = factor(rep(as.character(1:nrow(seqs)-1), ncol(seqs)), levels=as.character(1:nrow(seqs)-1))
df$ddG = factor(rep(rev(as.character(-1*(1:ncol(seqs)-1))), each=nrow(seqs)), levels=rev(as.character(-1*(1:ncol(seqs)-1))))

ggplot(df, aes_string(x="Mutations", y="ddG")) + 
  geom_tile(aes(fill=Score)) + 
  scale_fill_gradient(high="red", low="blue") + 
  geom_text(aes(label=round(Score, 2)), size = 2) + 
  ylab(expression(Delta*Delta*G)) + 
  theme_bw() + 
  theme(text=element_text(size=23, family="Helvetica"), axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0)), axis.line.x=element_line(color="black", size=1), 
        axis.ticks=element_line(color="black", size=1), panel.border=element_blank(),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1), axis.text.x.top = element_text(hjust=.5, size=15),
        axis.title.x.top = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("p53_heatmap_v1.pdf", width = 10, height=6)

# E-range to extract: -2, -4
dim(vals) = c(nrow(seqs), ncol(seqs))
p53.seqs = cbind(c(seqs[1, 13], seqs[2,11], seqs[3, 9]), c(vals[1, 13], vals[2,11], vals[3, 9]))

##################
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

# Cost function is input-output mapping
cost.fun = function(scored.list, base.scores) {
  return(sum(scored.list[[1]]))
}

# Start with optimal sequence
opt.seq = max.seq(NRLBModels, idx, 1)$BestSeq
GR.cost.fun.out = optimal.sequence(model.list = input.models, cost.function = cost.fun, seed.seq = opt.seq, verbose=TRUE)

plot.frontier(data = GR.cost.fun.out, derivative = F)
ggsave("GR_frontier_v1.pdf", width = 6, height=6)
plot.cost.range(data = GR.cost.fun.out) + xlab(expression(Delta*Delta*G)) + coord_flip()
ggsave("GR_costrange_v1.png", device = "png", width = 6, height=6)

seqs = GR.cost.fun.out$SequenceExamples
seqs[seqs=="NA"] = NA
seqs = seqs[,2:(ncol(seqs)-1)]
vals = score.seq.helper(sequences = seqs, base.seq = opt.seq, input.models)
vals = log(vals$Sum)
df = melt(vals, value.name="Score", varnames = c("Mutations", "Protein"))
df$Mutations = factor(rep(as.character(1:nrow(seqs)-1), ncol(seqs)), levels=as.character(1:nrow(seqs)-1))
df$ddG = factor(rep(rev(as.character(-.5*(1:ncol(seqs)-1))), each=nrow(seqs)), levels=rev(as.character(-.5*(1:ncol(seqs)-1))))

ggplot(df, aes_string(x="Mutations", y="ddG")) + 
  geom_tile(aes(fill=Score)) + 
  scale_fill_gradient(high="red", low="blue") + 
  geom_text(aes(label=round(Score, 2)), size = 2) + 
  ylab(expression(Delta*Delta*G)) + 
  theme_bw() + 
  theme(text=element_text(size=23, family="Helvetica"), axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0)), axis.line.x=element_line(color="black", size=1), 
        axis.ticks=element_line(color="black", size=1), panel.border=element_blank(),
        axis.text.x.bottom = element_text(angle = 45, hjust = 1), axis.text.x.top = element_text(hjust=.5, size=15),
        axis.title.x.top = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("GR_heatmap_v1.pdf", width = 10, height=6)


# E-range to extract: -2, -4
dim(vals) = c(nrow(seqs), ncol(seqs))
GR.seqs = cbind(c(seqs[1, 24], seqs[2,20], seqs[5, 16]), c(vals[1, 24], vals[2,20], vals[5, 16]))

save(p53.cost.fun.out, GR.cost.fun.out, file = "v1Run.rda")