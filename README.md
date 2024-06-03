# TP53 Reporter Optimization
This repository contains the code used to design and optimize the reporters in Trauernicht M, et al. (2023) Optimisation of TP53 reporters by systematic dissection of synthetic TP53 response elements. 	Nucleic Acids Research 51(18):9690-9702 (https://academic.oup.com/nar/article/51/18/9690/7256891).

[also contains some stuff for upcoming publication?]

[Include things related to ScoringTools?]

## Contents
* __p53_GR_initial_predictions.R__ Uses the RapidScore library and the associated sequence optimization code to build a set of sequences with varying affinity and number of mutations
* __response_element_scorer.R__ Uses the RapidScore library to compute affinities for the designed response elements
* __Images/__ contains plots made by `p53_GR_initial_predictions.R`
* __RData/__ contains the p53/GR models used to make predictions/score
* __data/__ contains files necessary for the two R scripts to run

mt20191119_TFs_list_updated.csv (Nov 2019):
Second follow-up list to Initial list of TFs considered to design reporters to study in mouse ES cells, updated list of TFs with corresponding motif and CIS-BP Motif ID

mt20191119_oligo_pool (3).fasta (Nov 2019):		
Draft .fasta file of the TF reporter oligo library for the 30 TFs in the updated TF list

mt20200619_oligo_pool_deep.fasta (Feb 2021): same sequences as in mt20210208_oligo_pool_deep_RE_only.fasta, but including the whole reporter (inclusive of the upstream adapter sequence, downstream minimal promoter, barcode and adapter)

mt20210208_oligo_pool_deep_RE_only.fasta (Feb 2021): contains all unique response elements (excluding primer adapters, minimal promoters, barcodes, etc.) for the designed P53 and GR reporter sequences used in the assay

fit.Zhang2017.GR.7297.json (Dec 2019):
ProBound GR fit 7297 used to create optimized sequences. This fit is derived from data generated as part of the Zhang 2017 paper.

* __out/__ contains files created by the R scripts
factor_info.xlsx (Nov 2019):	
File sent by Chaitanya to Max to help him design an initial set of promoters to explore (Nov 2019)	

p53scores/GRscores: files containing the outputs from running response_element_scorer.R

factor_info.xlsx (Nov 2019):	
File sent by Chaitanya to Max to help him design an initial set of promoters to explore (Nov 2019)	 

p53_GR_initial_predictions.RData:

## Requirements
R $`\geq `$ 4.2.1

The `RapidScore` library must be manually installed first - tarball is under `ScoringTools`

R packages:
- abind
- dplyr
- ggplot2
- reshape2
