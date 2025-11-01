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


# Installing  `R` Packages

The links below provide installation instructions and documentation for the `R` packages or code sets used to implement the six methods included in this repository:

* **`msm`**: See  [CRAN page](https://cran.r-project.org/web/packages/msm/index.html)

* **`msm-phase`**: Implemented within the `msm` package based on the method proposed by [Titman and Sharples (2010)](https://doi.org/10.1111/j.1541-0420.2009.01339.x)

* **`cthmm`**: See [R-Forge project page](https://r-forge.r-project.org/R/?group_id=1410) with manual available  [here](https://drive.google.com/file/d/1YapUT2xBFMzcjx9QxB_Cr3ejCfBLR4We/view?usp=sharing).  Example `R` code provided by the authors is available in the *Supporting Information* section of this [paper](https://onlinelibrary.wiley.com/doi/10.1111/biom.12252?msockid=03ab8c882dc66c933e759a8c2c146d2e)

* **`smms`**: See  [GitHub repository](https://github.com/NorskRegnesentral/smms)

* **`BayesTSM`**: See:  [GitHub repository](https://github.com/thomasklausch2/BayesTSM)

* **`hmm`**: The original code is available at [https://jonathanpw.github.io/research.html](https://jonathanpw.github.io/research.html)
  (Scroll to the *R code* for the 2020 published paper). However, the customized code included in this repository was provided directly by the authors.



