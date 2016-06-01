#' GEM: Fast association study for the interplay of Gene, Environment and Methylation
#'
#' The GEM package provides a highly efficient R tool suite for performing epigenome wide association
#' studies (EWAS). GEM provides three major functions named \code{\link{GEM_Emodel}},
#' \code{\link{GEM_Gmodel}} and \code{\link{GEM_GxEmodel}} to study the interplay of Gene,
#' Environment and Methylation (GEM). Within GEM, the existing "Matrix eQTL" package is
#' utilized and extended to study methylation quantitative trait loci (methQTL) and the
#' interaction of genotype and environment (GxE) to determine DNA methylation variation,
#' using matrix based iterative correlation and memory-efficient data analysis.
#' GEM can facilitate reliable genome-wide methQTL and GxE analysis on a standard laptop
#' computer within minutes.
#'
#' @examples
#'
#' ## Launch GEM GUI
#' #GEM_GUI()  # remove the hash symbol for running
#'
#' ## Checking the vignettes for more details
#' if(interactive()) browseVignettes(package = 'GEM')
#'
#' @seealso \code{\link{GEM_GUI}}
#' @references \url{https://github.com/fastGEM/GEM}
#' @author Hong Pan
#' @docType package
#' @name GEM-package
#'
NULL



#' Graphical User Interface (GUI) for GEM
#'
#' The user friendly GUI for runing GEM package easily and quickly
#'
#' The GEM package provides a highly efficient R tool suite for performing epigenome wide association
#' studies (EWAS). GEM provides three major functions named \code{\link{GEM_Emodel}},
#' \code{\link{GEM_Gmodel}} and \code{\link{GEM_GxEmodel}} to study the interplay of Gene,
#' Environment and Methylation (GEM). Within GEM, the pre-existing "Matrix eQTL" package is
#' utilized and extended to study methylation quantitative trait loci (methQTL) and the
#' interaction of genotype and environment (GxE) to determine DNA methylation variation,
#' using matrix based iterative correlation and memory-efficient data analysis.
#' GEM can facilitate reliable genome-wide methQTL and GxE analysis on a standard laptop
#' computer within minutes.
#'
#' @return GEM model analysis results
#' @import tcltk
#' @importFrom grDevices  dev.off jpeg
#' @importFrom graphics abline legend lines par plot points title
#' @importFrom methods new
#' @importFrom stats complete.cases lm pf pt qf qt
#' @importFrom utils flush.console tail tail.matrix write.table
#' @export
#' @seealso \code{\link{GEM-package}}
#' @examples
#' interactive()
#' #GEM_GUI()  ## remove the hash symbol to run
GEM_GUI <- function(){
    setModel <- TRUE   # model selection tag
    all_models <- c("Emodel", "Gmodel", "GxEmodel", "GplusEmodel", "GWASmodel")
    ifQuit <- FALSE

    while(setModel == TRUE){
        myModel <- modelSelection_GUI(all_models)

        if(myModel == "null"){
            setModel <- FALSE
            ifQuit <- TRUE
            break
        }else if(myModel == all_models[1]){
            paras <- modelPara_GUI(1, myModel, "Environment File Name", "Covariates File Name", "Methylation File Name")
            if(paras[[1]] == "no"){
                setModel <- FALSE
                if(paras[[8]] == "no"){
                    GEM_Emodel(env_file_name = paras[[3]],
                               covariate_file_name = paras[[4]],
                               methylation_file_name = paras[[5]],
                               Emodel_pv = as.numeric(paras[[6]]),
                               output_file_name = paste(paras[[2]], paste0(paras[[7]], "_results.txt"), sep = .Platform$file.sep),
                               qqplot_file_name = paste(paras[[2]], paste0(paras[[7]], "_qqplot.jpg"), sep = .Platform$file.sep))
                }else{
                    ifQuit <- TRUE
                }
            }

        }else if (myModel == all_models[2]){
            paras <- modelPara_GUI(5e-02, myModel, "SNP File Name", "Covariates File Name", "Methylation File Name")
            if(paras[[1]] == "no"){
                setModel <- FALSE
                if(paras[[8]] == "no"){
                    GEM_Gmodel(snp_file_name = paras[[3]],
                               covariate_file_name = paras[[4]],
                               methylation_file_name = paras[[5]],
                               Gmodel_pv = as.numeric(paras[[6]]),
                               output_file_name = paste(paras[[2]], paste0(paras[[7]], "_results.txt"), sep = .Platform$file.sep))
                }else{
                    ifQuit <- TRUE
                }
            }
        }else if(myModel == all_models[3]){
            paras <- modelPara_GUI(5e-02, myModel, "SNP File Name", "Covariates File Name", "Methylation File Name")
            if(paras[[1]] == "no"){
                setModel <- FALSE
                if(paras[[8]] == "no"){
                    GEM_GxEmodel(snp_file_name = paras[[3]],
                                 covariate_file_name = paras[[4]],
                                 methylation_file_name = paras[[5]],
                                 GxEmodel_pv = as.numeric(paras[[6]]),
                                 output_file_name = paste(paras[[2]], paste0(paras[[7]], "_results.txt"), sep = .Platform$file.sep))
                }else{
                    ifQuit <- TRUE
                }
            }
        }else if(myModel == all_models[4]){
            paras <- modelPara_GUI(5e-02, myModel, "SNP File Name", "Covariates File Name", "Methylation File Name")
            if(paras[[1]] == "no"){
                setModel <- FALSE
                if(paras[[8]] == "no"){
                    # GEM_GplusEmodel(snp_file_name = paras[[3]],
                    #                 covariate_file_name = paras[[4]],
                    #                 methylation_file_name = paras[[5]],
                    #                 GplusEmodel_pv = as.numeric(paras[[6]]),
                    #                 output_file_name = paste(paras[[2]], paste0(paras[[7]], "_results.txt"), sep = .Platform$file.sep),
                    #                 qqplot_file_name = paste(paras[[2]], paste0(paras[[7]], "_qqplot.jpg"), sep = .Platform$file.sep))
                }else{
                    ifQuit <- TRUE
                }
            }
        }else if(myModel == all_models[5]){
            paras <- modelPara_GUI(1, myModel, "Environment File Name", "SNP File Name", "Covariates File Name")
            if(paras[[1]] == "no"){
                setModel <- FALSE
                if(paras[[8]] == "no"){
                    GEM_GWASmodel(env_file_name = paras[[3]],
                                  snp_file_name = paras[[4]],
                                  covariate_file_name = paras[[5]],
                                  GWASmodel_pv = as.numeric(paras[[6]]),
                                  output_file_name = paste(paras[[2]], paste0(paras[[7]], "_results.txt"), sep = .Platform$file.sep),
                                  qqplot_file_name = paste(paras[[2]], paste0(paras[[7]], "_qqplot.jpg"), sep = .Platform$file.sep))
                }else{
                    ifQuit <- TRUE
                }
            }
        }
    }

    if(ifQuit){
        okMessage <- "Quit analysis."
    }else{
        okMessage <- paste0("Result file ",
                            paste0(paras[[7]], "_results.txt"),
                            " is under path:\n  ", paras[[2]])
        tkmessageBox(title = "GEM Package",
                     message = okMessage,
                     icon = "info", type = "ok")
    }

    message(okMessage)
}



modelPara_GUI <- function(pvalue, modelName = "Emodel",
                          fname1="env_file_name",
                          fname2="covariate_file_name",
                          fname3="methylation_file_name"){
    cur_dir <- getwd()
    changeM <- tclVar("no")     # model change label
    quitAlabel <- tclVar("no")
    dataDir <- tclVar(system.file('extdata',package='gem'))
    file1 <- tclVar("")
    file2 <- tclVar("")
    file3 <- tclVar("")
    pv <- tclVar(pvalue)
    projectName <- tclVar(paste0(modelName, "_analysis"))

    reset_dataDir <- function(){
        data_dir <- NULL
        data_dir <- tclvalue(tkchooseDirectory(title = "Choose your data dircetory ..."))
        if (!is.null(data_dir)) {
            tclvalue(dataDir) <- data_dir
            #tkmessageBox(title = "Data directory", message = paste0("Please make sure all your data is under path:\n  ", data_dir),icon = "info", type = "ok")
        }
    }

    reset_file1 <- function(){
        fn1 <- NULL
        fn1 <- tk_choose.files(default = paste(tclvalue(dataDir), "txt", sep = .Platform$file.sep),
                                  caption = "Select your file", multi = FALSE,
                                  filters = matrix(c("{txt files}", "{.txt}"),1, 2), index = 1)
        if (!is.null(fn1)) { tclvalue(file1) <- fn1 }
    }

    reset_file2 <- function(){
        fn2 <- NULL
        fn2 <- tk_choose.files(default = paste(tclvalue(dataDir), "txt", sep = .Platform$file.sep),
                                  caption = "Select your file", multi = FALSE,
                                  filters = matrix(c("{txt files}", "{.txt}"),1, 2), index = 1)
        if (!is.null(fn2)) { tclvalue(file2) <- fn2 }
    }

    reset_file3 <- function(){
        fn3 <- NULL
        fn3 <- tk_choose.files(default = paste(tclvalue(dataDir), "txt", sep = .Platform$file.sep),
                                  caption = "Select your file", multi = FALSE,
                                  filters = matrix(c("{txt files}", "{.txt}"),1, 2), index = 1)
        if (!is.null(fn3)) { tclvalue(file3) <- fn3 }
    }

    quitA <- function() {
        tclvalue(quitAlabel) <- "yes"
        tkdestroy(tt) }
    changeModelA <- function() {
        tclvalue(changeM) <- "yes"
        tkdestroy(tt)
    }

    submitA <- function() {
        has_error = FALSE
        if (as.numeric(tclvalue(pv)) > 1) {
            tkmessageBox(title = "GEM Package",
                         message = "Please provide correct cutoff P value (which must be less than 1).",
                         icon = "info", type = "ok")
            has_error = TRUE
        }

        if (has_error == FALSE) {
            tkdestroy(tt)
        }
    }

    box_length <- 50
    cell_width <- 3
    bt_width <- 8
    hb_width <- 2

    tt <- tktoplevel(borderwidth = 20)
    tkwm.title(tt, paste0(modelName, " Analysis"))

    projectName_label <- tklabel(tt, text = "Project Name :")
    projectName_entry <- tkentry(tt, textvariable = projectName, width = box_length)

    dataDir_label <- tklabel(tt, text = "Data Directory :")
    dataDir_entry <- tkentry(tt, textvariable = dataDir, width = box_length)
    dataDir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, command = reset_dataDir)

    file1_label <- tklabel(tt, text = paste0(fname1, " :"))
    file1_entry <- tkentry(tt, textvariable = file1, width = box_length)
    file1_button <- tkbutton(tt, text = " Select... ", width = bt_width, command = reset_file1)

    file2_label <- tklabel(tt, text = paste0(fname2, " :"))
    file2_entry <- tkentry(tt, textvariable = file2, width = box_length)
    file2_button <- tkbutton(tt, text = " Select... ", width = bt_width, command = reset_file2)

    file3_label <- tklabel(tt, text = paste0(fname3, " :"))
    file3_entry <- tkentry(tt, textvariable = file3, width = box_length)
    file3_button <- tkbutton(tt, text = " Select... ", width = bt_width, command = reset_file3)

    pv_label <- tklabel(tt, text = "Cutofff P value :")
    pv_entry <- tkentry(tt, textvariable = pv, width = box_length)

    quit_button <- tkbutton(tt, text = "Quit", command = quitA)
    changeModel_button <- tkbutton(tt, text = "Change model", command = changeModelA)
    submit_button <- tkbutton(tt, text = "Submit", command = submitA)


    tkgrid(projectName_label, projectName_entry, padx = cell_width)
    tkgrid.configure(projectName_label, projectName_entry ,sticky = "e")

    tkgrid(dataDir_label, dataDir_entry, dataDir_button, padx = cell_width)
    tkgrid.configure(dataDir_label, dataDir_entry, dataDir_button, sticky = "e")

    tkgrid(file1_label, file1_entry, file1_button, padx = cell_width)
    tkgrid.configure(file1_label, file1_entry, file1_button, sticky = "e")

    tkgrid(file2_label, file2_entry, file2_button, padx = cell_width)
    tkgrid.configure(file2_label, file2_entry, file2_button, sticky = "e")

    tkgrid(file3_label, file3_entry, file3_button, padx = cell_width)
    tkgrid.configure(file3_label, file3_entry, file3_button, sticky = "e")

    if(pvalue < 1){
        tkgrid(pv_label, pv_entry, padx = cell_width)
        tkgrid.configure(pv_label, pv_entry, sticky = "e")
    }

    tkgrid(tklabel(tt, text = "\n"), padx = cell_width)

    tkgrid(quit_button, changeModel_button, submit_button, padx = cell_width)
    tkgrid.configure(quit_button, sticky = "e")
    tkgrid.configure(submit_button, sticky = "w")

    tkwait.window(tt)

    return(list(tclvalue(changeM),
                tclvalue(dataDir),
                tclvalue(file1),
                tclvalue(file2),
                tclvalue(file3),
                tclvalue(pv),
                tclvalue(projectName),
                tclvalue(quitAlabel)
    ))

}


modelSelection_GUI <- function(models = c("Emodel", "Gmodel", "GxEmodel", "GplusEmodel", "GWASmodel")){

    all_models <- models
    myModel <- tclVar(all_models[1])
    cell_width <- 3

    quitM <- function(){
        tclvalue(myModel) <- "null"
        tkdestroy(tt)
    }
    helpM <- function() {
        tkmessageBox(title = "Model Information", message = "Emodel is to find the association between methylome and phenotype/environmental factors. The basic function is lm(methylation ~ env + covariates). The output is significance testing on environment.
                     \n\nGmodel is to find methyQTL that is methylation quantitative trait loci determined by genotyping. The basic function is lm(methylation ~ genotype + covariates). The output is significance testing on genotype.
                     \n\nGxEmodel is to explore how methylome is associated with the interplay between genotype and environmental factor.  The basic function is lm(methylation ~ genotype x environment +  covariates). The output is significance testing on genotype x environment."
                     ,icon = "info", type = "ok")}
    okM <- function() { tkdestroy(tt) }

    ## GUI
    tt <- tktoplevel(borderwidth = 20)
    tkwm.title(tt, "methylation_file_nameodel: GEM model selection")

    model_label <- tklabel(tt, text = "Please Select Your Analysis Model:")

    model_rbuts <- tkframe(tt)
    tkpack(tklabel(model_rbuts, text = "  "), side = "left")
    tkpack(tkradiobutton(model_rbuts, text = all_models[1], variable = myModel, value = all_models[1]), side = "left")
    tkpack(tkradiobutton(model_rbuts, text = all_models[2], variable = myModel, value = all_models[2]), side = "left")
    tkpack(tkradiobutton(model_rbuts, text = all_models[3], variable = myModel, value = all_models[3]), side = "left")
    #tkpack(tkradiobutton(model_rbuts, text = all_models[4], variable = myModel, value = all_models[4]), side = "left")
    #tkpack(tkradiobutton(model_rbuts, text = all_models[5], variable = myModel, value = all_models[5]), side = "left")

    lastline_buts <- tkframe(tt)
    tkpack(tkbutton(lastline_buts, text = "Quit", command = quitM), side = "left")
    tkpack(tklabel(lastline_buts, text = "          "), side = "left")
    tkpack(tkbutton(lastline_buts, text = "Help", command = helpM), side = "left")
    tkpack(tklabel(lastline_buts, text = "          "), side = "left")
    tkpack(tkbutton(lastline_buts, text = "OK", command = okM), side = "left")

    tkgrid(model_label, padx = cell_width)
    tkgrid.configure(model_label, sticky = "w")
    tkgrid(tklabel(tt, text = "\n"), padx = cell_width)

    tkgrid(model_rbuts, padx = cell_width)
    tkgrid.configure(model_rbuts, sticky = "e")
    tkgrid(tklabel(tt, text = "\n"), padx = cell_width)

    tkgrid(lastline_buts, padx = cell_width)

    tkwait.window(tt)

    return(tclvalue(myModel))
}


