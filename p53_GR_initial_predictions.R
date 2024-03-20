######################
### Load libraries ###
######################
source("ScoringTools/base/WLDOS_SequenceOptimizer.R")
source("ScoringTools/base/PlotLibrary.R")

####################
### p53 Analysis ###
####################

# Load NRLB p53 model and some attributes
load('RData/p53_models.RData')

# Cost function is the identity
cost.fun = function(scored.list, base.scores) {
  return(scored.list[[1]][1])
}

# Run Optimizations - Start with optimal sequence
p53.cost.fun.out = optimal.sequence(model.list = p53.model, cost.function = cost.fun, seed.seq = opt.p53.seq, verbose=TRUE)

# Visualize
plot.frontier(data = p53.cost.fun.out, derivative = F)
ggsave("p53_frontier_v1.pdf", width = 6, height=6)
plot.cost.range(data = p53.cost.fun.out) + xlab(expression(Delta*Delta*G)) + coord_flip()
ggsave("Images/p53_costrange_v1.png", device = "png", width = 6, height=6)

seqs = p53.cost.fun.out$SequenceExamples
seqs[seqs=="NA"] = NA
vals = score.seq.helper(sequences = seqs, base.seq = opt.p53.seq, p53.model)
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
ggsave("Images/p53_heatmap_v1.pdf", width = 10, height=6)

# E-range to extract: -2, -4
dim(vals) = c(nrow(seqs), ncol(seqs))
p53.seqs = data.frame(Seq = c(seqs[1, 13], seqs[4,11], seqs[3, 9]), 
                      ddG = c(vals[1, 13], vals[4,11], vals[3, 9]), 
                      Affinity=exp(c(vals[1, 13], vals[4,11], vals[3, 9])))

###################
### GR Analysis ###
###################
# Load GR and get betas from a probound model (7297)
load('RData/GR_models.RData')

GR.cost.fun.out = optimal.sequence(model.list = GR.model, cost.function = cost.fun, seed.seq = opt.GR.seq, verbose=TRUE)

#Visualize
plot.frontier(data = GR.cost.fun.out, derivative = F)
ggsave("GR_frontier_v1.pdf", width = 6, height=6)
plot.cost.range(data = GR.cost.fun.out) + xlab(expression(Delta*Delta*G)) + coord_flip()
ggsave("Images/GR_costrange_v1.png", device = "png", width = 6, height=6)

seqs = GR.cost.fun.out$SequenceExamples
seqs[seqs=="NA"] = NA
seqs = seqs[,2:(ncol(seqs)-1)]
vals = score.seq.helper(sequences = seqs, base.seq = opt.GR.seq, GR.model)
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
ggsave("Images/GR_heatmap_v1.pdf", width = 10, height=6)


# E-range to extract: -2, -4
dim(vals) = c(nrow(seqs), ncol(seqs))
GR.seqs = data.frame(Seq = c(seqs[1, 24], seqs[2,20], seqs[5, 16]), 
                     ddG = c(vals[1, 24], vals[2,20], vals[5, 16]), 
                     Affinity=exp(c(vals[1, 24], vals[2,20], vals[5, 16])))

# Store actual condensed sequence outs into RData file
save(p53.cost.fun.out, GR.cost.fun.out, p53.seqs, GR.seqs, file = "RData/p53_GR_initial_predictions.RData")
