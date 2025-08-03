# Modelling Skewed and Heavy-Tailed Errors in Bayesian Mediation Analysis

<p align="center">
  <img src="/img/left.png" alt="Conditional total effect — plot 1" width="45%">
  <img src="/img/right.png" alt="Conditional total effect — plot 2" width="45%">
</p>


This repository contains all material needed to reproduce the results reported in **“Modelling Skewed and Heavy-Tailed Errors in Bayesian Mediation Analysis”** (Li, Steel, & Zhang, 2025).

The methods described in the article are implemented in the R package **_FlexBayesMed_**. The package is under active development; new features—such as one-line estimation of complex moderated-mediation models—are added regularly. Future updates will appear here. Feel free to contact the maintainer with questions or suggestions.

## Repository layout

Here is a brief description of the repository structure.

| Folder          | Contents                                                                    |
| --------------- | ----------------------------------------------------------------------------|
| `FlexBayesMed/` | Source code for the R package _FlexBayesMed_.                               |
| `papercode/`    | Scripts for reproducing analyses and figures in Li, Steel, & Zhang (2025). |
| `document/`     | PDF reference manual generated from the package documentation.              |
| `demo/`         | User guide (PDF) plus synthetic example data sets.                          |
| `data/`         | Datasets used in the case study in the paper and in the demo.               |


## Key features of **_FlexBayesMed_**

1. **Distribution utilities**  
   * Compute the PDF, CDF, quantiles, and Fisher's moment coefficients of skewness of the *centred two-piece Student-t* distribution.  
   * Draw random samples from the same distribution.

2. **Flexible linear regression** — [`flex_mr()`]  
   Fit linear models with skewed or heavy-tailed errors using any of the following error distributions:  
   * centred two-piece Student-t  
   * skew-normal  
   * Student-t  
   * normal

3. **Bayesian single-level mediation analysis** — [`flex_med()`]  
   * Handles the same four error distributions as `flex_mr()`.  
   * Performs Bayes-factor testing of the mediation effect.  
   * Lets you customise the hyperparameters of priors for $\gamma$ (skewness parameter) and $\nu$ (tail parameter).

4. **Coming soon**  
   One-line wrappers for a range of moderated-mediation models and other extensions are in active development.

> **Tip**  
> You can already fit complex models—such as moderated mediation—by applying `flex_mr()` to each regression component separately. See the worked example in `demo/`.

## Dependencies

**_FlexBayesMed_** relies on three core packages:

| Package | Purpose |
|---------|---------|
| **rstan** | Fits models with HMC/NUTS |
| **bridgesampling** | Computes marginal likelihoods & Bayes factors based on the Stan objects|
| **Rmpfr** | Provides high-precision arithmetic |


## How to cite
