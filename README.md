# BayesianNetworkCoextinction
This repository contains the data, code, figures, and tables for "Maximum trophic level predicts food webs’ susceptibility to coextinctions" (Li et al.)

Authors: \*Henry Li, Anna Ekl&ouml;f, &dagger;Gy&ouml;rgy Barab&aacute;s, &dagger;Laura E. Dee

&dagger; denotes co-senior author

\*Corresponding Author: Henry Li, henry.li-1@colorado.edu

This repository is archived in Zenodo at https://doi.org/10.5281/zenodo.21939116

The preprint of the associated manuscript is available in bioRxiv at

## Replication Instructions
Set working directory to your local repository location and run the following `R` scripts in the order listed.

1. `code/generate_webs.R`, which creates the following data objects:
   - `data/stratified_webs.rds`
   - `data/niche_webs.rds`
   - `data/empirical_webs.rds`
   - `data/empirical_webs_ES.rds`
2. `code/run_simulations.R`, which creates the following data objects:
   - `data/stratified_outcome.rds` and `data/stratified_outcome_small.rds`
   - `data/niche_outcome.rds` and `data/niche_outcome_small.rds`
   - `data/empirical_outcome.rds` and `data/empirical_outcome_small.rds`
   - `data/empirical_ES_outcome.rds` and `data/empirical_ES_outcome_small.rds`
3. `code/fit_parameter.R`, which creates the following data objects:
   - `data/stratified_param_fit.rds`
   - `data/niche_param_fit.rds` and `data/niche_param_adjusted.rds`
   - `data/empirical_param_fit.rds`
   - `data/empirical_omnivory_deviation.rds`
4. `code/visualizations.R`, which creates figures found in `figures_and_tables/`

## Figure and Table Key
Descriptions for figures and tables found in `figures_and_tables/`

1. Main text:
   - Figure 1: Comparison of extinction scenario simulation approaches (a and b), and food web augmented with ecosystem service (ES) nodes and links (c)
   - Figure 2: Excess Extinction and Fragility as a function of direct extinction risk in stratified food webs
   - Figure 3: Excess Extinction as a function of direct extinction risk in niche model food webs with different connectance values
   - Figure 4: Comparison of true maximum trophic level and effective maximum trophic level in niche model and empirical food webs, and extinction simulation outcomes for empirical food webs compared to corresponding model food webs
   - Figure 5: Excess Loss of Ecosystem Services, grouped by empirical food web
   - Table 1: Description of calculated quantities
2. Supporting Information:
   - Figure S1: Comparison of food web models
   - Figure S2: Excess Extinction and Fragility as a function of food web size in niche model food webs
   - Figure S3: Excess Extinction and Fragility as a function of food web connectance in niche model food webs
   - Figure S4: Excess Extinction and Fragility as a function of food web size and connectance in empirical food webs
   - Figure S5: Excess Extinction and Fragility as a function of food web size in empirical food webs
   - Figure S6: Excess Extinction and Fragility as a function of food web connectance in empirical food webs
   - Figure S7: Maximum trophic level in niche model food webs as a function of food web structure
   - Figure S8: Omnivory index of empirical food webs
   - Figure S9: Comparison of maximum trophic level in empirical food webs and corresponding niche model food webs
   - Figure S10: Excess Extinction as a function of maximum trophic level in stratified food webs
   - Figure S11: Fragility as a function of maximum trophic level in stratified food webs
   - Figure S12: Excess Extinction as a function of maximum trophic level in niche model food webs
   - Figure S13: Fragility as a function of maximum trophic level in niche model food webs
   - Figure S14: Excess Extinction as a function of maximum trophic level in empirical food webs
   - Figure S15: Fragility as a function of maximum trophic level in empirical food webs
   - Figure S16: Extinction simulation outcomes for empirical food webs compared to corresponding model food webs
   - Table S1: Description of empirical food webs used for analysis
