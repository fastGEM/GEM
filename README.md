GEM
----------------------
Gene, Environment and Methylation (GEM): A tool suite to efficiently navigate large scale epigenome wide association studies and integrate genotype and interaction between genotype and environment

## Installation

### Install from Bioconductor

Currently, the GEM package has been submitted to Bioconductor and under review. Once the package is on Bioconductor, you can install the package with following commands:

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEM")
```

### Install from github

```R
## install dependent packages
if(!require("methods")){
    install.packages("methods")
}
if(!require("ggplot2")){
    install.packages("ggplot2")
}
if(!require("Rcpp")){
    install.packages("Rcpp")
}
if(!require("digest")){
    install.packages("digest")
}
if(!require("devtools")){
    install.packages("devtools")
}

## install gem package
devtools::install_github("fastGEM/GEM")
```

## Launch the GUI

```R
library(GEM)
GEM_GUI()
```


## Example codes

```R
library(GEM)
DATADIR = system.file('extdata',package='GEM')
RESULTDIR = getwd()

env_file_name = paste(DATADIR, "env.txt", sep = .Platform$file.sep)
covariates_file_name = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
covariates_file_name_gxe = paste(DATADIR, "gxe.txt", sep = .Platform$file.sep)
methylation_file_name = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
snp_file_name = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)

Emodel_pv = 1
Gmodel_pv = 1e-04
GxEmodel_pv = 1

Emodel_result_file_name = paste(RESULTDIR, "Result_Emodel.txt", sep = .Platform$file.sep)
Gmodel_result_file_name = paste(RESULTDIR, "Result_Gmodel.txt", sep = .Platform$file.sep)
GxEmodel_result_file_name = paste(RESULTDIR, "Result_GxEmodel.txt", sep = .Platform$file.sep)

Emodel_qqplot_file_name = paste(RESULTDIR, "QQplot_Emodel.jpg", sep = .Platform$file.sep)

GEM_Emodel(env_file_name, covariates_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name,Emodel_qqplot_file_name)

GEM_Gmodel(snp_file_name, covariates_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name)

GEM_GxEmodel(snp_file_name, covariates_file_name_gxe, methylation_file_name, GxEmodel_pv, GxEmodel_result_file_name)

```
