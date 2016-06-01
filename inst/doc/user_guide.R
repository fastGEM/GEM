## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(dpi = 300)
knitr::opts_chunk$set(cache=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("GEM")

## ---- results='hide'-----------------------------------------------------
require(GEM)

## ---- eval=FALSE---------------------------------------------------------
#  GEM_GUI()

## ------------------------------------------------------------------------
DATADIR = system.file('extdata',package='GEM')
dir(path = DATADIR)

## ----echo=FALSE----------------------------------------------------------
read.table(paste(DATADIR, "cov.txt", sep = .Platform$file.sep), header = TRUE, row.names = 1)[1, 1:10]

## ----echo=FALSE----------------------------------------------------------
read.table(paste(DATADIR, "env.txt", sep = .Platform$file.sep), header = TRUE, row.names = 1)[1, 1:10]

## ----echo=FALSE----------------------------------------------------------
read.table(paste(DATADIR, "gxe.txt", sep = .Platform$file.sep), header = TRUE, row.names = 1)[1:2, 1:10]

## ---- echo = FALSE-------------------------------------------------------
read.table(paste(DATADIR, "methylation.txt", sep = .Platform$file.sep), header = TRUE)[1:5, 1:7]

## ---- echo = FALSE-------------------------------------------------------
read.table(paste(DATADIR, "snp.txt", sep = .Platform$file.sep), header = TRUE)[1:5, 1:15]

## ---- eval=TRUE, fig.width=6, fig.height=6, out.width=600, out.height=600----
env_file_name = paste(DATADIR, "env.txt", sep = .Platform$file.sep)
covariate_file_name = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
methylation_file_name = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
Emodel_pv = 1
Emodel_result_file_name = "Result_Emodel.txt"
Emodel_qqplot_file_name = "QQplot_Emodel.jpg"
GEM_Emodel(env_file_name, covariate_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name, Emodel_qqplot_file_name, savePlot=FALSE)

## ----echo=FALSE----------------------------------------------------------
head(read.table(paste(getwd(), "Result_Emodel.txt", sep = .Platform$file.sep), header = TRUE))

## ---- eval=TRUE----------------------------------------------------------
snp_file_name = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)
covariate_file_name = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
methylation_file_name = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
Gmodel_pv = 1e-04
Gmodel_result_file_name = "Result_Gmodel.txt"

GEM_Gmodel(snp_file_name, covariate_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name)

## ----echo=FALSE----------------------------------------------------------
head(read.table(paste(getwd(), "Result_Gmodel.txt", sep = .Platform$file.sep), header = TRUE))

## ---- eval=TRUE, fig.keep='all', dev='png', fig.width=8, fig.height=5, out.width=600, out.height=400----
snp_file_name = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)
covariate_file_name = paste(DATADIR, "gxe.txt", sep = .Platform$file.sep)
methylation_file_name = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
GxEmodel_pv = 1
GxEmodel_result_file_name = "Result_GxEmodel.txt"
GEM_GxEmodel(snp_file_name, covariate_file_name, methylation_file_name, GxEmodel_pv, GxEmodel_result_file_name, topKplot = 1, savePlot=FALSE)

## ----echo=FALSE----------------------------------------------------------
head(read.table(paste(getwd(), "Result_GxEmodel.txt", sep = .Platform$file.sep), header = TRUE))

