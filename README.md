## Modelling the implementation of narrow versus broader spectrum antibiotics in the empiric treatment of E. coli bacteraemia

This repo contains code to reproduce the analyses and figures for "Modelling the implementation of narrow versus broader spectrum antibiotics in the empiric treatment of E. coli bacteraemia".

## Environment setup:
R version: 4.2.2

The following libraries are used:
```R
library(cowplot)
library(deSolve)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)
library('lhs')
library('sensitivity')
```


## To reproduce the results:
1. Run ```Baseline_analysis.R``` to get the results from the baseline analysis. For the baseline analysis, to explore uncertainty, a normal distribution was fit to several highly uncertain or setting-dependent parameters (i.e., treatment duration [T], empiric therapy duration [1/δ], breakthrough resistance rates [⍺], the effect of inappropriate therapy [𝜀], population-level effects of antibiotic use on resistance levels [], and baseline mortality rates [D]), which we denote variable parameters. The mean for each variable parameter was set to the sourced baseline parameter value, with standard deviations (SD) set such that the 95% CI of the distribution was within the ranges used for the sensitivity analyses, 
2. Run ```Sensitivity_analysis.R``` to get results from the sensitivity analyses. For sensitivity analyses, we split variable parameter values into low, medium, and high values.
3. Run ```Multivariate_LHS_analysis_analysis.R``` to get Latin hypercube sampling results. Multivariate sensitivity analyses were conducted using Latin Hypercube Sampling (1000 samples). Partial rank correlation coefficients with 95% CIs were calculated for each parameter.


## Sharing/Access information
SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC0 1.0 Universal (CC0 1.0) Public Domain

2. Was data derived from another source? NA


## Paper Citation
TBD
