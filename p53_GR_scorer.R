library(NRLBtools)
library(jsonlite)
library(RapidScore)
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
  # Load seqs (switch to appropriate input file)
#  seqs = read.table("/Users/chaitanya/Documents/Research/SequenceEngineering/NKI/mt20210208_oligo_pool_deep_RE_only.fasta", 
#                    sep=",", header=F, stringsAsFactors = F)
  seqs = read.table("/Users/chaitanya/Documents/Research/SequenceEngineering/NKI/mt20200619_oligo_pool_deep.fasta", 
                    sep=",", header=F, stringsAsFactors = F)
  names(seqs) = c("Name", "Seq")
  # Load new sequences
  load("~/Documents/Research/SequenceEngineering/NKI/mt20210308_hg38_promoters-500bp.RData")
  load("~/Documents/Research/SequenceEngineering/NKI/mt20210412_p53_REs.RData")
}


# Load and score p53
idx = which(NRLBModelInfo$Protein=="p53" & NRLBModelInfo$Info=="WildType")
energy.adjustment = log(max.seq(NRLBModels, idx, 1)$MaxAffinity)
p53.scorer = new(RapidScore, -30)
p53.scorer$add(NRLBModels[[2]][[idx]][[1]]$NB, NRLBModels[[2]][[idx]][[1]]$DB)
p53.seqs = dplyr::filter(seqs, grepl("Trp53", Name))
p53.seqs = cbind(p53.seqs, SumAffinity=sapply(p53.seqs$Seq, FUN=function(x) sum(exp(p53.scorer$scoreBulkChar(x)[[1]]-energy.adjustment))))
write.table(x = p53.seqs, file = "~/Desktop/p53scores.tsv", sep = "\t", col.names = T, quote = F, row.names = F)

# Create scored dfs
p53.seqs = dplyr::filter(seqs, grepl("Trp53", Name))
per.pos = sapply(p53.seqs$Seq, FUN = function(x) exp(p53.scorer$scoreBulkChar(x)[[1]]-energy.adjustment))
names(per.pos) = p53.seqs$Name
per.pos = list(Info=p53.seqs, Scores=per.pos)
save(per.pos,file = "~/Documents/Research/SequenceEngineering/NKI/p53_per_position.Rda")

# Score genomic tss
per.pos = vector(mode="list", length=nrow(tss_annotated))
for (i in 1:length(per.pos)) {
  per.pos[[i]] = exp(p53.scorer$scoreBulkChar(tss_annotated$seq[i])[[1]]-energy.adjustment)
}
names(per.pos) = tss_annotated$transcript_id
tss.per.pos = list(Info=tss_annotated, Scores=per.pos)
save(tss.per.pos, file = "~/Documents/Research/SequenceEngineering/NKI/p53_tss_per_position.Rda")

# Score response elements
per.pos = vector(mode="list", length=nrow(p53_REs_export))
for (i in 1:length(per.pos)) {
  per.pos[[i]] = exp(p53.scorer$scoreBulkChar(p53_REs_export$seq[i])[[1]]-energy.adjustment)
}
names(per.pos) = p53_REs_export$name
re.per.pos = list(Info=p53_REs_export, Scores=per.pos)
save(re.per.pos, file = "~/Documents/Research/SequenceEngineering/NKI/p53_REs_per_position.Rda")

# Load and score GR
json.con = file(description = "fit.Zhang2017.GR.7297.json", open="r")
GR.json = jsonlite::stream_in(json.con)
close(json.con)
# Extract betas
GR.bm = GR.json$coefficients$bindingModes[[1]]
mono = GR.bm$mononucleotide[[2]]
di = as.numeric(GR.bm$dinucleotide[[2]])
gr.scorer = new(RapidScore, -30)
gr.scorer$add(mono, di)
gr.seqs = dplyr::filter(seqs, grepl("Gr", Name))
gr.seqs = cbind(gr.seqs, SumAffinity=sapply(gr.seqs$Seq, FUN=function(x) sum(exp(gr.scorer$scoreBulkChar(x)[[1]]))))
write.table(x = gr.seqs, file = "~/Desktop/grscores.tsv", sep = "\t", col.names = T, quote = F, row.names = F)
