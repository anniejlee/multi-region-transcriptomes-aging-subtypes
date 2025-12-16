# Multi-Region Transcriptomic Subtyping 
This repository contains R code for the manuscript entitled:  
**“Multi-region brain transcriptomes uncover two subtypes of aging individuals with differences in Alzheimer’s disease risk and the impact of APOE ε4”.**

This repository provides simplified R scripts illustrating the analytic workflow used in the manuscript.  
The pipeline includes:

- Sparse Multi-View Canonical Correlation Analysis (MultiCCA)  
- Region-specific K-means clustering  
- Cross-region meta-clustering using Non-negative Matrix Factorization (NMF)  

All examples are based on toy simulated datasets and do not include real biological data.

---------------------------------------------------------------
## Software Requirements
---------------------------------------------------------------
These scripts require only standard R packages:

- **PMA** (for MultiCCA)  
- **NMF** (for NMF meta-clustering)

Install in R with:

```r
install.packages("PMA")
install.packages("NMF")
```

---------------------------------------------------------------
## Description of R Scripts
---------------------------------------------------------------
- **dis-cluster.R** : This script performs the discovery cohort analysis using toy data. It runs Sparse MultiCCA, region-specific K-means clustering, and NMF meta-clustering.
- **rep-cluster.R** : This script performs the replication cohort analysis. It projects replication samples into the discovery CCA space, assigns cluster labels, and applies NMF meta-clustering.


