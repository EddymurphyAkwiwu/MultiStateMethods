# Repository Overview

This repository contains several folders for statistical methods for modelling multi-state cancer progression described in the paper: *"A comparison of methods for modelling multi-state cancer progression using screening data with censoring after intervention."*. These folders are as follows: 

- **msm**
- **msm-phase**
- **cthmm**
- **smms**
- **BayesTSM**
- **hmm**
- **Proof_of_CDF**

## Description of Folders

### **1. Folders for R Software Package Implementations**
The first six folders (**msm**, **msm-phase**, **cthmm**, **smms**, **BayesTSM**, and **hmm**) each contain sample `R` code implementing the corresponding R packages described in the above paper:  


Each of these folders has two subfolders: **Exponential** and **Weibull**. These subfolders include:

1. Simulated data based on a 3-state model, assuming the true cumulative distribution function (CDF) for progression times generated in a large population setting (n = 10^6). The data are based on two settings: p = 2 (i.e., two covariates) and strong censoring mechanism, and p = 2 and medium censoring mechanism (see the above paper for details).
2. `R` code for generating data (e.g., n = 1000 or n = 2000) based on the two settings mentioned above.
3. `R` code implementing the method for the respective `R` package.

Each of the above methods is implemented under the setting with p = 2 and strong censoring mechanism, with the p = 2 and medium censoring setting commented out in the `R` file.  To use, please uncomment.


For other simulation settings, see the above paper for parameter settings.  Additionally, the **controllers** and **source** subfolders, located within each subfolder, contain `R` functions essential for generating the simulated data.

### **2. Proof_of_CDF Folder**
The **Proof_of_CDF** folder contains R code to demonstrate how the cumulative distribution function (CDF) (also known as the cumulative incidence function, CIF) can be derived using functions from the **msm** (and **msm-phase**) package. Specifically, it demonstrates the use of the following functions:

- `ppass.msm()`
- `qmatrix.msm()`
- `pmatrix.msm()`

As an example, the demonstration assumes an exponentially distributed 3-state model without covariates (i.e., p = 0). For other models, such as Weibull distributions or those with covariates, a similar approach applies. Please refer to the **msm** package user manual for guidance on fitting models with covariates.
