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
* __out/__ contains files created by the R scripts

## Requirements
R $`\geq `$ 4.2.1

The `RapidScore` library must be manually installed first - tarball is under `ScoringTools`

R packages:
- abind
- dplyr
- ggplot2
- reshape2
