# TP53 Reporter Optimization
This repository contains the code used to design and optimize the reporters in Trauernicht M, et al. (2023) Optimisation of TP53 reporters by systematic dissection of synthetic TP53 response elements. 	Nucleic Acids Research 51(18):9690-9702 (https://academic.oup.com/nar/article/51/18/9690/7256891).

[also contains some stuff for upcoming publication?]

## Contents
* `p53_GR_initial_predictions.R` Uses the RapidScore library and the associated sequence optimization code to build a set of sequences with varying affinity and number of mutations
* `response_element_scorer.R` Uses the RapidScore library to compute affinities for the designed response elements
* `Images/` contains plots made by `p53_GR_initial_predictions.R`
* `RData/` contains the p53/GR models used to make predictions/score
* `data/` contains input files necessary for the two R scripts to run
  * `mt20191119_TFs_list_updated.csv` Second follow-up list to Initial list of TFs considered to design reporters to study in mouse ES cells, updated list of TFs with corresponding motif and CIS-BP Motif ID
  * `mt20210208_oligo_pool_deep_RE_only.fasta` contains all unique response elements (excluding primer adapters, minimal promoters, barcodes, etc.) for the designed P53 and GR reporter sequences used in the assay
  * `mt20200619_oligo_pool_deep.fasta` same sequences as in mt20210208_oligo_pool_deep_RE_only.fasta, but including the whole reporter (inclusive of the upstream adapter sequence, downstream minimal promoter, barcode and adapter)
  * `fit.Zhang2017.GR.7297.json` [ProBound](https://www.nature.com/articles/s41587-022-01307-0) GR fit 7297 used to create optimized sequences. This fit is derived from data generated as part of the Zhang 2017 paper.
* `out/` contains files created by the R scripts
  * `factor_info.xlsx` contains manually-modified motifs from `mt20191119_TFs_list_updated.csv` that are designed to eliminate binding for the given transcription factor. __For new publication__
  * `GR/p53scores_RE/_full.tsv` files output by `response_element_scorer.R` that contain the GR/p53 affinity scores for sequences in `mt20210208_oligo_pool_deep_RE_only.fasta`/`mt20200619_oligo_pool_deep.fasta` 
  * `p53_GR_initial_predictions.RData` RData files containing the outputs of the sequence optimization runs
  * `FullSequences_v1.xslx` contains the optimized sequences running to create an affinity ladder for p53/GR.

## Requirements
R $`\geq `$ 4.2.1

The `RapidScore` library must be manually installed first - tarball is under `ScoringTools`

R packages:
- abind
- dplyr
- ggplot2
- reshape2
