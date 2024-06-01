# TP53 Reporter Optimization
This repository contains the code used to design and optimize the reporters in Trauernicht M, et al. (2023) Optimisation of TP53 reporters by systematic dissection of synthetic TP53 response elements. 	Nucleic Acids Research 51(18):9690-9702 (https://academic.oup.com/nar/article/51/18/9690/7256891).

## Contents
p53_GR_initial_predictions.R: R file that uses the RapidScore library and the associated sequence optimization code to build a set of sequences with varying affinity and number of mutations

response_element_scorer.R: R file that uses the RapidScore library to compute affinities for the designed response elements

Images/: stores plots made by p53_GR_initial_predictions.R

RData/: stores the dataframes and other output of p53_GR_initial_predictions.R and the p53/GR models used to make them

data/: stores the info used to build the models

out/: stores the output from code

## Requirements
R $`\geq `$ 4.2.1

The `RapidScore` library must be manually installed first - tarball is under `ScoringTools`

R packages:
abind
dplyr
ggplot2
reshape2
