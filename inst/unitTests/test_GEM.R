test_GEM <- function() {
  DATADIR = system.file('extdata',package='GEM')
  RESULTDIR = getwd()
  
  env_file_name = paste(DATADIR, "env.txt", sep = .Platform$file.sep)
  covariates_file_name = paste(DATADIR, "cov.txt", sep = .Platform$file.sep)
  covariates_file_name_gxe = paste(DATADIR, "gxe.txt", sep = .Platform$file.sep)
  methylation_file_name = paste(DATADIR, "methylation.txt", sep = .Platform$file.sep)
  snp_file_name = paste(DATADIR, "snp.txt", sep = .Platform$file.sep)
  
  Emodel_pv = 1
  Gmodel_pv = 1e-04
  GxEmodel_pv = 1e-4
  
  Emodel_result_file_name = paste(RESULTDIR, "Result_Emodel.txt", sep = .Platform$file.sep)
  Gmodel_result_file_name = paste(RESULTDIR, "Result_Gmodel.txt", sep = .Platform$file.sep)
  GxEmodel_result_file_name = paste(RESULTDIR, "Result_GxEmodel.txt", sep = .Platform$file.sep)
  
  Emodel_qqplot_file_name = paste(RESULTDIR, "QQplot_Emodel.jpg", sep = .Platform$file.sep)
  
  GEM_Emodel(env_file_name, covariates_file_name, methylation_file_name, Emodel_pv, Emodel_result_file_name,
             Emodel_qqplot_file_name)
  GEM_Gmodel(snp_file_name, covariates_file_name, methylation_file_name, Gmodel_pv, Gmodel_result_file_name)
  GEM_GxEmodel(snp_file_name, covariates_file_name_gxe, methylation_file_name, GxEmodel_pv,
               GxEmodel_result_file_name)
  
  checkTrue(length(dir(path = RESULTDIR)) > 1)
}

