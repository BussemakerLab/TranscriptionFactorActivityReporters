#############################
### Load libraries & data ###
#############################
source("ScoringTools/base/WLDOS_SequenceOptimizer.R")
source("ScoringTools/base/PlotLibrary.R")
# Load NRLB p53 model and some attributes
load('RData/p53_models.RData')
# Load GR and get betas from a probound model (7297)
load('RData/GR_models.RData')
# Load seqences 
seqs.RE = read.table("data/mt20210208_oligo_pool_deep_RE_only.fasta", sep=",", header=F, stringsAsFactors = F)
seqs.full = read.table("data/mt20200619_oligo_pool_deep.fasta", sep=",", header=F, stringsAsFactors = F)
names(seqs.RE) = names(seqs.full) = c("Name", "Seq")

####################
### p53 Analysis ###
####################
# Use rapid score to compute affinities
p53.scorer = new(RapidScore, -30)
p53.scorer$add(p53.model$p53$NB, p53.model$p53$DB)

# Score response elements and full sequences
seqs = dplyr::filter(seqs.RE, grepl("Trp53", Name))
scores = cbind(seqs, SumAffinity=sapply(seqs$Seq, FUN=function(x) sum(exp(p53.scorer$scoreBulkChar(x)[[1]]))))
write.table(x = scores, file = "out/p53scores_RE.tsv", sep = "\t", col.names = T, quote = F, row.names = F)

seqs = dplyr::filter(seqs.full, grepl("Trp53", Name))
scores = cbind(seqs, SumAffinity=sapply(seqs$Seq, FUN=function(x) sum(exp(p53.scorer$scoreBulkChar(x)[[1]]))))
write.table(x = scores, file = "out/p53scores_full.tsv", sep = "\t", col.names = T, quote = F, row.names = F)

###################
### GR Analysis ###
###################
# Use rapid score to compute affinities
GR.scorer = new(RapidScore, -30)
GR.scorer$add(GR.model$GR$NB, GR.model$GR$DB)

# Score response elements and full sequences
seqs = dplyr::filter(seqs.RE, grepl("Gr", Name))
scores = cbind(seqs, SumAffinity=sapply(seqs$Seq, FUN=function(x) sum(exp(p53.scorer$scoreBulkChar(x)[[1]]))))
write.table(x = scores, file = "out/GRscores_RE.tsv", sep = "\t", col.names = T, quote = F, row.names = F)

seqs = dplyr::filter(seqs.full, grepl("Gr", Name))
scores = cbind(seqs, SumAffinity=sapply(seqs$Seq, FUN=function(x) sum(exp(p53.scorer$scoreBulkChar(x)[[1]]))))
write.table(x = scores, file = "out/GRscores_full.tsv", sep = "\t", col.names = T, quote = F, row.names = F)
