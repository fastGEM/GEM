GEM
----------------------


## Installation

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

## install gem package
install.packages("devtools")
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

env_file_name = paste(DATADIR, "env.txt", sep = "")
covariates_file_name = paste(DATADIR, "cov.txt", sep = "")
covariates_file_name_gxe = paste(DATADIR, "gxe.txt", sep = "")
methylation_file_name = paste(DATADIR, "methylation.txt", sep = "")
snp_file_name = paste(DATADIR, "snp.txt", sep = "")

Emodel_pv = 1
Gmodel_pv = 1e-04
GxEmodel_pv = 1

Emodel_result_file_name = paste(RESULTDIR, "Result_Emodel.txt", sep = "")
Gmodel_result_file_name = paste(RESULTDIR, "Result_Gmodel.txt", sep = "")
GxEmodel_result_file_name = paste(RESULTDIR, "Result_GxEmodel.txt", sep = "")

Emodel_qqplot_file_name = paste(RESULTDIR, "QQplot_Emodel.jpg", sep = "")
Gmodel_qqplot_file_name = paste(RESULTDIR, "QQplot_Gmodel.jpg", sep = "")
GxEmodel_qqplot_file_name = paste(RESULTDIR, "QQplot_GxEmodel.jpg", sep = "")

GEM_Emodel(env_file_name, covariates_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name,
    Emodel_qqplot_file_name)
GEM_Gmodel(snp_file_name, covariates_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name,
    Gmodel_qqplot_file_name)
GEM_GxEmodel(snp_file_name, covariates_file_name_gxe, methylation_file_name, GxEmodel_pv,
    GxEmodel_result_file_name)

```
