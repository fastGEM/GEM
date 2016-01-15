##################### DATA

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

LM_Emodel_result_file_name = paste(RESULTDIR, "LM_Result_Emodel.txt", sep = "")
LM_Gmodel_result_file_name = paste(RESULTDIR, "LM_Result_Gmodel.txt", sep = "")
LM_GxEmodel_result_file_name = paste(RESULTDIR, "LM_Result_GxEmodel.txt", sep = "")

########################## Benchmarking
env <- read.table(env_file_name, header = T, sep = "\t", check.names = F)
envdata <- matrix(unlist(env), ncol = dim(env)[2])
envdata <- envdata[, 2:dim(envdata)[2]]

cov <- read.table(covariates_file_name, header = T, sep = "\t", check.names = F)
covdata <- matrix(unlist(cov), ncol = dim(cov)[2])
len <- dim(cov)[2] - 1
x <- rep(1, len)
y <- cov[, 2:dim(cov)[2]]
covdata <- t(rbind(x, y))

methydata <- read.table(methylation_file_name, header = T, sep = "\t")
snpdata <- read.table(snp_file_name, header = T, sep = "\t")
numCpG <- dim(methydata)[1]
numSNP <- dim(snpdata)[1]

# emodel
start1 = proc.time()[3]
myresult <- matrix(rep(NA, numCpG * 4), ncol = 4)
for (i in 1:numCpG) {
    cpg1 <- unlist(methydata[i, 2:dim(methydata)[2]])
    fit1 <- lm(cpg1 ~ envdata + covdata[, ])
    f1 <- summary(fit1)
    pv <- f1$coefficient[2, 4]
    if (pv <= Emodel_pv) {
        tmp1 <- round(c(f1$coefficient[2, 1], f1$coefficient[2, 2], -log10(f1$coefficient[2, 
            4])), 6)
        myresult[i, ] <- cbind(as.character(methydata[i, 1]), t(tmp1))
    }
}
colnames(myresult) <- c("CpGID", "Est", "SD", "neglogpv")
write.table(myresult, LM_Emodel_result_file_name, row.name = F, sep = "\t", quote = F)
end1 = proc.time()[3]
cat("Time for Emodel: ", end1 - start1, " seconds\n")


#### gmodel
start = proc.time()[3]
for (i in 1:numCpG) {
    start1 = proc.time()[3]
    cpg1 <- unlist(methydata[i, 2:dim(methydata)[2]])
    for (j in 1:numSNP) {
        snp1 <- unlist(snpdata[j, 2:2:dim(methydata)[2]])
        fit1 <- lm(cpg1 ~ snp1 + covdata[, 1:dim(covdata)[2]])
        f1 <- summary(fit1)
        (pv <- f1$coefficient[2, 4])
        if (pv <= Gmodel_pv) {
            tmp1 <- round(c(f1$coefficient[2, 1], f1$coefficient[2, 2], -log10(f1$coefficient[2, 
                4])), 6)
            myresult1 <- cbind(as.character(methydata[i, 1]), as.character(snpdata[j, 
                1]), t(tmp1))
            colnames(myresult1) <- c("CpGID", "snpid", "Est", "SD", "neglogpv")
            if (i == 1) 
                write.table(myresult1, LM_Gmodel_result_file_name, row.name = F, sep = "\t", 
                  quote = F) else {
                write.table(myresult1, LM_Gmodel_result_file_name, row.name = F, col.name = F, 
                  sep = "\t", append = T, quote = F)
            }
        }
    }
    end1 = proc.time()[3]
    # cat ('finished ', i, ' at ', end1-start1, ' seconds\n');
}
end = proc.time()[3]
cat("time for Gmodel:", end - start, "seconds\n")




# gxemodel
start = proc.time()[3]
for (i in 1:numCpG) {
    cpg1 <- unlist(methydata[i, 2:dim(methydata)[2]])
    for (j in 1:numSNP) {
        snp1 <- unlist(snpdata[j, 2:2:dim(methydata)[2]])
        fit1 <- lm(cpg1 ~ snp1 * envdata + covdata)
        f1 <- summary(fit1)
        if (dim(f1$coefficient)[1] == 5) {
            tmp1 <- round(c(f1$coefficient[5, 1], f1$coefficient[5, 2], -log10(f1$coefficient[5, 
                4])), 6)
            myresult1 <- cbind(as.character(methydata[i, 1]), as.character(snpdata[j, 
                1]), t(tmp1))
            colnames(myresult1) <- c("CpGID", "snpid", "B12Est", "B12SD", "B12neglogpv")
            if (i == 1) 
                write.table(myresult1, LM_GxEmodel_result_file_name, row.name = F, 
                  sep = "\t") else {
                write.table(myresult1, LM_GxEmodel_result_file_name, row.name = F, 
                  col.name = F, sep = "\t", append = T)
            }
        }
        
    }
}
end = proc.time()[3]
cat("time for GxEmodel:", end - start, "seconds\n")
 
