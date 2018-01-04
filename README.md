# Riemannian distances and divergences

## Introduction
This code accompanies a the publication "Assessment of Riemannian distances and divergences for SSVEP-based BCI" which aims at assessing the impact of several distances/divergences on a real EEG dataset.

## Important files
- `main.m`

  Classifies SSVEP trials from 12 subjects, using MDM with distances/divergences described in the aforementioned publication: arithmetic, harmonic, Riemannian, log-euclid, Kullback-Leibler, S-divergence, $\alpha$-divergence, Bhattacharrya, Wasserstein, and Jeffreys.
  The main output of this file are the classification accuracies and the computation time for each method/metric

- `alpha_cross_validation.m`  

  In the $\alpha$-divergence, the value of $\alpha$ is deternined through cross-validation.

- `loaddata.m`  

  Implement a function called in `main.m`. The path to the dataset is hardcoded herein.

- `swelling_effect_analysis.m`

  Computes the determinants and traces of means of covariance matrices computed with different distances/divergences.
  It also computes them for each individual covariance matrix used in the computation of the means.  

  Only 1 subject and one class are used for illustrative purposes.  
