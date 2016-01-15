#' Emodel Analysis
#'
#' Emodel is to find the association between methylome and phenotype/environmental factors. The basic function is lm(methylation ~ env + covariates). The output is significance testing on environment.
#'
#' @param env_file_name Phenotype/environmental factor data file
#' @param covariates_file_name Covariates data file
#' @param methylation_file_name Methylation data file
#' @param Emodel_pv Pvalue-cutoff (default is 1.0)
#' @param output_file_name output file name
#' @param qqplot_file_name qqplot file name
#'
#' @return save results automatically
#' @export
#'
#' @examples
#' DATADIR = system.file('extdata',package='gem')
#' RESULTDIR = getwd()
#' env_file_name = paste(DATADIR, "env.txt", sep = "")
#' covariates_file_name = paste(DATADIR, "cov.txt", sep = "")
#' covariates_file_name_gxe = paste(DATADIR, "gxe.txt", sep = "")
#' methylation_file_name = paste(DATADIR, "methylation.txt", sep = "")
#' snp_file_name = paste(DATADIR, "snp.txt", sep = "")
#' Emodel_pv = 1
#' Emodel_result_file_name = paste(RESULTDIR, "Result_Emodel.txt", sep = "")
#' Emodel_qqplot_file_name = paste(RESULTDIR, "QQplot_Emodel.jpg", sep = "")
#' GEM_Emodel(env_file_name, covariates_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name, Emodel_qqplot_file_name)
GEM_Emodel <-
    function(env_file_name, covariates_file_name, methylation_file_name,
             Emodel_pv, output_file_name, qqplot_file_name) {

        errorCovariance = numeric();

        env <- SlicedData$new();
        env$fileDelimiter = "\t";      # the TAB character
        env$fileOmitCharacters = "NA"; # denote missing values;
        env$fileSkipRows = 1;          # one row of column labels
        env$fileSkipColumns = 1;       # one column of row labels
        env$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        env$LoadFile(env_file_name);

        cvrt <- SlicedData$new();
        cvrt$fileDelimiter = "\t";      # the TAB character
        cvrt$fileOmitCharacters = "NA"; # denote missing values;
        cvrt$fileSkipRows = 1;          # one row of column labels
        cvrt$fileSkipColumns = 1;       # one column of row labels
        if (length(covariates_file_name) > 0) {
            cvrt$LoadFile(covariates_file_name);
        }

        cpg = SlicedData$new();
        cpg$fileDelimiter = "\t";      # the TAB character
        cpg$fileOmitCharacters = "NA"; # denote missing values;
        cpg$fileSkipRows = 1;          # one row of column labels
        cpg$fileSkipColumns = 1;       # one column of row labels
        cpg$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        cpg$LoadFile(methylation_file_name);

        ## Run the analysis
        Emodel <- Matrix_eQTL_engine2(
            snps = env,
            gene = cpg,
            cvrt = cvrt,
            output_file_name =  NULL,
            pvOutputThreshold = Emodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE,
            addInfo = "CpGs"
        )

        unlink(output_file_name);
        ## Results:
        cat('Analysis done in: ', Emodel$time.in.sec, ' seconds', '\n');
        #show(Emodel$all$eqtls)
        #R2 = Emodel$all$eqtls$statistic ^ 2 / (Emodel$all$eqtls$statistic ^ 2 + Emodel$param$dfFull);
        result_Emodel <- cbind(
            as.character(Emodel$all$eqtls$gene),
            Emodel$all$eqtls$beta,
            Emodel$all$eqtls$statistic,
            round(-log10(Emodel$all$eqtls$pvalue),6)
        )
        colnames(result_Emodel) <- c("cpg", "beta", "stats", "logpv")
        write.table(
            result_Emodel, output_file_name, sep = "\t", row.names = F, quote = F
        )

        jpeg(qqplot_file_name, width = 2000, height = 2000, res = 300)
        plot(Emodel, pch = 16, cex = 0.7);
        dev.off()
    }


#' Gmodel Analysis
#'
#' Gmodel is to find methyQTL that is methylation quantitative trait loci determined by genotyping. The basic function is lm(methylation ~ genotype + covariates). The output is significance testing on genotype.
#'
#' @param snp_file_name Genotype data file
#' @param covariates_file_name Covariates data file
#' @param methylation_file_name Methylation data file
#' @param Gmodel_pv Pvalue-cutoff(default is 1.0e-5 in most cases)
#' @param output_file_name output file name
#' @param qqplot_file_name qqplot file name
#'
#' @return save results automatically
#' @export
#'
#' @examples
#' DATADIR = system.file('extdata',package='gem')
#' RESULTDIR = getwd()
#' env_file_name = paste(DATADIR, "env.txt", sep = "")
#' covariates_file_name = paste(DATADIR, "cov.txt", sep = "")
#' covariates_file_name_gxe = paste(DATADIR, "gxe.txt", sep = "")
#' methylation_file_name = paste(DATADIR, "methylation.txt", sep = "")
#' snp_file_name = paste(DATADIR, "snp.txt", sep = "")
#' Gmodel_pv = 1e-04
#' Gmodel_result_file_name = paste(RESULTDIR, "Result_Gmodel.txt", sep = "")
#' Gmodel_qqplot_file_name = paste(RESULTDIR, "QQplot_Gmodel.jpg", sep = "")
#' GEM_Gmodel(snp_file_name, covariates_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name, Gmodel_qqplot_file_name)
GEM_Gmodel <-
    function(snp_file_name, covariates_file_name, methylation_file_name,
             Gmodel_pv, output_file_name)
    {

        errorCovariance = numeric();

        snp = SlicedData$new();
        snp$fileDelimiter = "\t";      # the TAB character
        snp$fileOmitCharacters = "NA"; # denote missing values;
        snp$fileSkipRows = 1;          # one row of column labels
        snp$fileSkipColumns = 1;       # one column of row labels
        snp$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        snp$LoadFile(snp_file_name);

        cvrt = SlicedData$new();
        cvrt$fileDelimiter = "\t";      # the TAB character
        cvrt$fileOmitCharacters = "NA"; # denote missing values;
        cvrt$fileSkipRows = 1;          # one row of column labels
        cvrt$fileSkipColumns = 1;       # one column of row labels
        if (length(covariates_file_name) > 0) {
            cvrt$LoadFile(covariates_file_name);
        }

        cpg = SlicedData$new();
        cpg$fileDelimiter = "\t";      # the TAB character
        cpg$fileOmitCharacters = "NA"; # denote missing values;
        cpg$fileSkipRows = 1;          # one row of column labels
        cpg$fileSkipColumns = 1;       # one column of row labels
        cpg$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        cpg$LoadFile(methylation_file_name);


        ## Run the analysis
        Gmodel = Matrix_eQTL_engine2(
            snps = snp,
            gene = cpg,
            cvrt = cvrt,
            output_file_name = NULL,
            pvOutputThreshold = Gmodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE,
            addInfo = "methQTL"
        );

        unlink(output_file_name);
        ## Results:
        cat('Analysis done in: ', Gmodel$time.in.sec, ' seconds', '\n');
        R2 = Gmodel$all$eqtls$statistic ^ 2 / (Gmodel$all$eqtls$statistic ^ 2 + Gmodel$param$dfFull);
        result_Gmodel <- cbind(
            as.character(Gmodel$all$eqtls$gene),
            as.character(Gmodel$all$eqtls$snps),
            Gmodel$all$eqtls$beta,
            Gmodel$all$eqtls$statistic,
            round(-log10(Gmodel$all$eqtls$pvalue),6)
        )
        colnames(result_Gmodel) <- c("cpg", "snp", "beta", "stats", "logpv")
        write.table(
            result_Gmodel, output_file_name, sep = "\t", row.names = F, quote = F
        )
    }




#' GxEmodel
#'
#' GxEmodel is to explore how methylome is associated with the interplay between genotype and environmental factor.  The basic function is lm(methylation ~ genotype x environment +  covariates). The output is significance testing on genotype x environment.
#'
#' @param snp_file_name Genotype data file
#' @param covariates_file_name Covariates data file
#' @param methylation_file_name Methylation data file
#' @param GxEmodel_pv Pvalue-cutoff(default is 1.0e-5 in most cases)
#' @param output_file_name output file name
#' @param topKplot plot the top k number of assocations, if NULL, not plot will be saved
#'
#' @return save results automatically
#' @import ggplot2
#' @export
#'
#' @examples
#' DATADIR = system.file('extdata',package='gem')
#' RESULTDIR = getwd()
#' env_file_name = paste(DATADIR, "env.txt", sep = "")
#' covariates_file_name = paste(DATADIR, "cov.txt", sep = "")
#' covariates_file_name_gxe = paste(DATADIR, "gxe.txt", sep = "")
#' methylation_file_name = paste(DATADIR, "methylation.txt", sep = "")
#' snp_file_name = paste(DATADIR, "snp.txt", sep = "")
#' GxEmodel_pv = 1
#' GxEmodel_result_file_name = paste(RESULTDIR, "Result_GxEmodel.txt", sep = "")
#' GxEmodel_qqplot_file_name = paste(RESULTDIR, "QQplot_GxEmodel.jpg", sep = "")
#' GEM_GxEmodel(snp_file_name, covariates_file_name_gxe, methylation_file_name, GxEmodel_pv, GxEmodel_result_file_name)
GEM_GxEmodel <-
    function(snp_file_name, covariates_file_name, methylation_file_name,
             GxEmodel_pv, output_file_name, topKplot = 10)
    {

        errorCovariance = numeric();

        snp = SlicedData$new();
        snp$fileDelimiter = "\t";      # the TAB character
        snp$fileOmitCharacters = "NA"; # denote missing values;
        snp$fileSkipRows = 1;          # one row of column labels
        snp$fileSkipColumns = 1;       # one column of row labels
        snp$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        snp$LoadFile(snp_file_name);

        cvrt = SlicedData$new();
        cvrt$fileDelimiter = "\t";      # the TAB character
        cvrt$fileOmitCharacters = "NA"; # denote missing values;
        cvrt$fileSkipRows = 1;          # one row of column labels
        cvrt$fileSkipColumns = 1;       # one column of row labels
        if (length(covariates_file_name) > 0) {
            cvrt$LoadFile(covariates_file_name);
        }

        cpg = SlicedData$new();
        cpg$fileDelimiter = "\t";      # the TAB character
        cpg$fileOmitCharacters = "NA"; # denote missing values;
        cpg$fileSkipRows = 1;          # one row of column labels
        cpg$fileSkipColumns = 1;       # one column of row labels
        cpg$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        cpg$LoadFile(methylation_file_name);


        ## Run the analysis
        GxEmodel = Matrix_eQTL_engine2(
            snps = snp,
            gene = cpg,
            cvrt = cvrt,
            output_file_name = NULL,
            pvOutputThreshold = GxEmodel_pv,
            useModel = modelLINEAR_CROSS,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE,
            addInfo = "cpg-snp pairs"
        );

        ## Results:
        cat('Analysis done in: ', GxEmodel$time.in.sec, ' seconds', '\n');
        R2 = GxEmodel$all$eqtls$statistic ^ 2 / (GxEmodel$all$eqtls$statistic ^
                                                     2 + GxEmodel$param$dfFull);

        if(!is.null(topKplot) && is.numeric(topKplot) && topKplot > 0){
            if(topKplot > length(GxEmodel$all$eqtl$pvalue))
                topKplot <- length(GxEmodel$all$eqtl$pvalue)
            for(i in 1:topKplot){
                snpiname <- GxEmodel$all$eqtl$snps[i]
                cpginame <- GxEmodel$all$eqtl$gene[i]

                snpi<- snp$FindRow(snpiname)$row[1, ]
                cpgi<- cpg$FindRow(cpginame)$row[1, ]

                sceData <- data.frame(snp = factor(snpi, levels=c(1,2,3), labels=c("AA", "AB", "BB")),
                                      cpg = as.numeric(cpgi),
                                      env = as.numeric(cvrt$getSlice(1)[dim(cvrt)[1], ]))

                dp <- ggplot(sceData, aes(x=env, y=cpg)) + geom_point() +
                    geom_smooth(method=lm, se=FALSE, colour = "red") +
                    facet_grid(. ~ snp) + theme_bw() + xlab("Environmental factor") + ylab(paste0(cpginame, " profile")) +
                    theme(strip.text.x = element_text(size=12, face="bold"))
                ggsave(paste0(dirname(output_file_name), .Platform$file.sep,
                              snpiname, " x ", cpginame, ".pdf"),
                       dp, width=12, height=7)
            }
        }

        result_GxEmodel <- cbind (
            as.character(GxEmodel$all$eqtls$gene),
            as.character(GxEmodel$all$eqtls$snps),
            GxEmodel$all$eqtls$beta,
            GxEmodel$all$eqtls$statistic,
            round(-log10(GxEmodel$all$eqtls$pvalue),6)
            )

        colnames(result_GxEmodel) <- c("cpg", "snp", "beta", "stats", "logpv")
        write.table(
            result_GxEmodel, output_file_name, sep = "\t", row.names = F, quote = F
        )

        #jpeg(qqplot_file_name, width = 2000, height = 2000, res = 300)
        #plot(GxEmodel, pch = 16, cex = 0.7);
        #dev.off()
}




#' GplusEmodel
#'
#' @param snp_file_name
#' @param covariates_file_name
#' @param methylation_file_name
#' @param GplusEmodel_pv
#' @param output_file_name
#' @param qqplot_file_name
#'
#' @return save results automatically
#' @export
GEM_GplusEmodel <-
    function(snp_file_name, covariates_file_name, methylation_file_name,
             GplusEmodel_pv, output_file_name, qqplot_file_name)
    {

        errorCovariance = numeric();

        snp = SlicedData$new();
        snp$fileDelimiter = "\t";      # the TAB character
        snp$fileOmitCharacters = "NA"; # denote missing values;
        snp$fileSkipRows = 1;          # one row of column labels
        snp$fileSkipColumns = 1;       # one column of row labels
        snp$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        snp$LoadFile(snp_file_name);

        cvrt = SlicedData$new();
        cvrt$fileDelimiter = "\t";      # the TAB character
        cvrt$fileOmitCharacters = "NA"; # denote missing values;
        cvrt$fileSkipRows = 1;          # one row of column labels
        cvrt$fileSkipColumns = 1;       # one column of row labels
        if (length(covariates_file_name) > 0) {
            cvrt$LoadFile(covariates_file_name);
        }

        cpg = SlicedData$new();
        cpg$fileDelimiter = "\t";      # the TAB character
        cpg$fileOmitCharacters = "NA"; # denote missing values;
        cpg$fileSkipRows = 1;          # one row of column labels
        cpg$fileSkipColumns = 1;       # one column of row labels
        cpg$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        cpg$LoadFile(methylation_file_name);


        ## Run the analysis
        GplusEmodel = Matrix_eQTL_engine2(
            snps = snp,
            gene = cpg,
            cvrt = cvrt,
            output_file_name = NULL,
            pvOutputThreshold = GplusEmodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE
        );

        ## Results:
        cat('Analysis done in: ', GplusEmodel$time.in.sec, ' seconds', '\n');
        R2 = GplusEmodel$all$eqtls$statistic ^ 2 / (GplusEmodel$all$eqtls$statistic ^
                                                        2 + GplusEmodel$param$dfFull);
        result_GplusEmodel <- cbind(
            as.character(GplusEmodel$all$eqtl$gene),
            as.character(GplusEmodel$all$eqtl$snps),
            round(-log10(GplusEmodel$all$eqtl$pvalue),6),
            round(GplusEmodel$all$eqtl$FDR, 6),
            round(R2,6)
        )
        colnames(result_GplusEmodel) <- c("cpg", "snp", "logpv", "FDR", "R2")
        write.table(
            result_GplusEmodel, output_file_name, sep = "\t", row.names = F, quote =
                F
        )

        jpeg(
            qqplot_file_name, width = 2000, height = 2000, res = 300
        )
        plot(GplusEmodel, pch = 16, cex = 0.7);
        dev.off()
    }



#' GWASmodel
#'
#' @param env_file_name
#' @param snp_file_name
#' @param covariates_file_name
#' @param GWASmodel_pv
#' @param output_file_name
#' @param qqplot_file_name
#'
#' @return save results automatically
#' @export
GEM_GWASmodel <-
    function(env_file_name, snp_file_name, covariates_file_name,
             GWASmodel_pv, output_file_name, qqplot_file_name)
    {

        errorCovariance = numeric();

        snp = SlicedData$new();
        snp$fileDelimiter = "\t";      # the TAB character
        snp$fileOmitCharacters = "NA"; # denote missing values;
        snp$fileSkipRows = 1;          # one row of column labels
        snp$fileSkipColumns = 1;       # one column of row labels
        snp$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        snp$LoadFile(snp_file_name);

        cvrt = SlicedData$new();
        cvrt$fileDelimiter = "\t";      # the TAB character
        cvrt$fileOmitCharacters = "NA"; # denote missing values;
        cvrt$fileSkipRows = 1;          # one row of column labels
        cvrt$fileSkipColumns = 1;       # one column of row labels
        if (length(covariates_file_name) > 0) {
            cvrt$LoadFile(covariates_file_name);
        }

        env = SlicedData$new();
        env$fileDelimiter = "\t";      # the TAB character
        env$fileOmitCharacters = "NA"; # denote missing values;
        env$fileSkipRows = 1;          # one row of column labels
        env$fileSkipColumns = 1;       # one column of row labels
        env$fileSliceSize = 2000;      # read file in slices of 2,000 rows
        env$LoadFile(env_file_name);


        ## Run the analysis
        GWASmodel = Matrix_eQTL_engine2(
            snps = snp,
            gene = env,
            cvrt = cvrt,
            output_file_name = NULL,
            pvOutputThreshold = GWASmodel_pv,
            useModel = modelLINEAR,
            errorCovariance = errorCovariance,
            verbose = FALSE,
            pvalue.hist = "qqplot",
            min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = FALSE,
            addInfo = "SNPs"
        );

        ## Results:
        cat('Analysis done in: ', GWASmodel$time.in.sec, ' seconds', '\n');
        R2 = GWASmodel$all$eqtls$statistic ^ 2 / (GWASmodel$all$eqtls$statistic ^
                                                      2 + GWASmodel$param$dfFull);
        result_GWASmodel <- cbind(
            as.character(GWASmodel$all$eqtl$snps),
            round(-log10(GWASmodel$all$eqtl$pvalue),6),
            round(GWASmodel$all$eqtl$FDR, 6),
            round(R2,6)
        )
        colnames(result_GWASmodel) <- c("snp", "logpv", "FDR", "R2")
        write.table(
            result_GWASmodel, output_file_name, sep = "\t", row.names = F, quote = F
        )

        jpeg(
            qqplot_file_name, width = 2000, height = 2000, res = 300
        )
        plot(GWASmodel, pch = 16, cex = 0.7);
        dev.off()
    }


#' Emodel
#'
#' @param Y
#' @param X
#' @param covariates
#'
#' @return save results automatically
#' @export
Emodel <- function(Y,X,covariates) {
    #Y=bX+cov
    U1 = crossprod(covariates, Y)
    U2 = solve(crossprod(covariates), U1)
    ytr = Y - covariates %*% U2

    U3 = crossprod(covariates, X)
    U4 = solve(crossprod(covariates), U3)
    Xtr = X - covariates %*% U4

    k = dim(covariates)[2]
    n = dim(covariates)[1]

    df = n - k - 2
    b = as.vector(crossprod(ytr, Xtr) / colSums(Xtr ^ 2))
    Xtr2 = colSums(Xtr ^ 2)
    sig = (sum(ytr ^ 2) - b ^ 2 * Xtr2) / df
    err = sqrt(sig * (1 / Xtr2))
    #p = 2 * pnorm(-abs(b/err))
    p = 2 * pt(-abs(b / err), df)
    logp = -log10(p)

    result <- cbind(round(b, 6), round(err, 6), round(logp,6))
    return(result)
}
