##################### DATA
source("GEM_model.R")

DATADIR = "data/"
RESULTDIR = "RESULT/"

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

quit()

########################## Benchmarking
env <- read.table(env_file_name, header = T, sep = "\t", check.names = F)
envdata <- matrix(unlist(env), ncol = dim(env)[2])
envdata <- envdata[, 2:dim(envdata)[2]]

cov <- read.table(covariates_file_name_gxe, header = T, sep = "\t", check.names = F)
covdata <- matrix(unlist(cov), ncol = dim(cov)[2])
len <- dim(cov)[2] - 1
x <- rep(1, len)
y <- cov[, 2:dim(cov)[2]]
covdata <- t(rbind(x, y))

methydata <- read.table(methylation_file_name, header = T, sep = "\t")
snpdata <- read.table(snp_file_name, header = T, sep = "\t")
numCpG <- dim(methydata)[1]

LM_Emodel_result_file_name = paste(RESULTDIR, "VitD_nmol_L_LM_Emodel_result.txt", sep = "")
LM_Gmodel_result_file_name = paste(RESULTDIR, "VitD_nmol_L_LM_Gmodel_result.txt", sep = "")
LM_GxEmodel_result_file_name = paste(RESULTDIR, "VitD_nmol_L_LM_GxEmodel_result.txt",
    sep = "")
LM_GplusEmodel_result_file_name = paste(RESULTDIR, "VitD_nmol_L_LM_GplusEmodel_result.txt",
    sep = "")
LM_GWASmodel_result_file_name = paste(RESULTDIR, "VitD_nmol_L_LM_GWASmodel_result.txt",
    sep = "")

# emodel
start = proc.time()[3]
# numCpG=2
for (i in 1:numCpG) {
    cpg1 <- unlist(methydata[i, 2:dim(methydata)[2]])
    # snp1<- unlist(snpdata[9, 2:2:dim(methydata)[2] ]) because this snp does not have
    # missing value.
    fit <- lm(cpg1 ~ envdata + covdata[, 1:7])
    f <- summary(fit)
    tmp <- round(c(f$coefficient[2, 1], f$coefficient[2, 2], -log10(f$coefficient[2,
        4])), 6)
    myresult <- cbind(as.character(methydata[i, 1]), t(tmp))
    colnames(myresult) <- c("CpGID", "B12Est", "B12SD", "B12neglogpv")
    if (i == 1)
        write.table(myresult, LM_Emodel_result_file_name, row.name = F, sep = "\t",
            quote = F) else {
        write.table(myresult, LM_Emodel_result_file_name, row.name = F, col.name = F,
            sep = "\t", append = T, quote = F)
    }
}  #check the Emodel result in 'VitD_nmol_L_Emodel_result.txt' .
end = proc.time()[3]
cat("time for Emodel:", end - start, "seconds\n")

# gmodel
start = proc.time()[3]
numSNP <- dim(snpdata)[1]
# numSNP=2
for (i in 1:numCpG) {
    cpg1 <- unlist(methydata[i, 2:dim(methydata)[2]])
    for (j in 1:numSNP) {
        snp1 <- unlist(snpdata[j, 2:2:dim(methydata)[2]])
        fit1 <- lm(cpg1 ~ snp1 + covdata[, 1:7])
        f1 <- summary(fit1)
        tmp1 <- round(c(f1$coefficient[2, 1], f1$coefficient[2, 2], -log10(f1$coefficient[2,
            4])), 6)
        myresult1 <- cbind(as.character(methydata[i, 1]), as.character(snpdata[j, 1]),
            t(tmp1))
        colnames(myresult1) <- c("CpGID", "snpid", "B12Est", "B12SD", "B12neglogpv")
        if (i == 1)
            write.table(myresult1, LM_Gmodel_result_file_name, row.name = F, sep = "\t",
                quote = F) else {
            write.table(myresult1, LM_Gmodel_result_file_name, row.name = F, col.name = F,
                sep = "\t", append = T, quote = F)
        }
    }
}
end = proc.time()[3]
cat("time for Gmodel:", end - start, "seconds\n")
# check the Gmodel result in 'VitD_nmol_L_Gmodel_result.txt' .

# gxemodel
start = proc.time()[3]
for (i in 1:numCpG) {
    cpg1 <- unlist(methydata[i, 2:dim(methydata)[2]])
    for (j in 1:numSNP) {
        snp1 <- unlist(snpdata[j, 2:2:dim(methydata)[2]])
        fit1 <- lm(cpg1 ~ snp1 * envdata + covdata[, 1:7])
        f1 <- summary(fit1)
        tmp1 <- round(c(f1$coefficient[10, 1], f1$coefficient[10, 2], -log10(f1$coefficient[10,
            4])), 6)
        myresult1 <- cbind(as.character(methydata[i, 1]), as.character(snpdata[j, 1]),
            t(tmp1))
        colnames(myresult1) <- c("CpGID", "snpid", "B12Est", "B12SD", "B12neglogpv")
        if (i == 1)
            write.table(myresult1, LM_GxEmodel_result_file_name, row.name = F, sep = "\t") else {
            write.table(myresult1, LM_GxEmodel_result_file_name, row.name = F, col.name = F,
                sep = "\t", append = T)
        }
    }
}  #check the GxEmodel result in 'VitD_nmol_L_GxEmodel_result.txt' .
end = proc.time()[3]
cat("time for GxEmodel:", end - start, "seconds\n")

# gplusemodel??????I report G not E????
start = proc.time()[3]
for (i in 1:numCpG) {
    cpg1 <- unlist(methydata[i, 2:dim(methydata)[2]])
    for (j in 1:numSNP) {
        snp1 <- unlist(snpdata[j, 2:2:dim(methydata)[2]])
        fit1 <- lm(cpg1 ~ snp1 + covdata[, 1:8])
        f1 <- summary(fit1)
        tmp1 <- round(c(f1$coefficient[2, 1], f1$coefficient[2, 2], -log10(f1$coefficient[2,
            4])), 6)
        myresult1 <- cbind(as.character(methydata[i, 1]), as.character(snpdata[j, 1]),
            t(tmp1))
        colnames(myresult1) <- c("CpGID", "snpid", "B12Est", "B12SD", "B12neglogpv")
        if (i == 1)
            write.table(myresult1, LM_GplusEmodel_result_file_name, row.name = F, sep = "\t") else {
            write.table(myresult1, LM_GplusEmodel_result_file_name, row.name = F, col.name = F,
                sep = "\t", append = T)
        }
    }
}  #check the GplusEmodel result in 'VitD_nmol_L_GplusEmodel_result.txt' .
end = proc.time()[3]
cat("time for GplusEmodel:", end - start, "seconds\n")

start = proc.time()[3]
# GWASmodel
for (j in 1:numSNP) {
    snp1 <- unlist(snpdata[j, 2:2:dim(methydata)[2]])
    fit1 <- lm(envdata ~ snp1 + covdata[, c(1, 2, 6, 7)])
    f1 <- summary(fit1)
    tmp1 <- round(c(f1$coefficient[2, 1], f1$coefficient[2, 2], -log10(f1$coefficient[2,
        4])), 6)
    myresult1 <- cbind(as.character(snpdata[j, 1]), t(tmp1))
    colnames(myresult1) <- c("snpid", "B12Est", "B12SD", "B12neglogpv")
    # print(myresult1)
    if (j == 1)
        write.table(myresult1, LM_GWASmodel_result_file_name, row.name = F, sep = "\t") else {
        write.table(myresult1, LM_GWASmodel_result_file_name, row.name = F, col.name = F,
            sep = "\t", append = T)
    }

}
# check the GplusEmodel result in 'VitD_nmol_L_GWASmodel_result.txt' .
end = proc.time()[3]
cat("time for GWASmodel:", end - start, "seconds\n")





