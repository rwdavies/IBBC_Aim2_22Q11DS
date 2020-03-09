nget_and_sanitize <- function(what) {
   gsub("\\\\", "", Sys.getenv(what))
}

R_dir <- get_and_sanitize("R_DIR")
results_dir <- get_and_sanitize("RESULTS_DIR")
external_dir <- get_and_sanitize("EXTERNAL_DIR")
iBBC <- get_and_sanitize("IBBC")
pheno_file_name_no_extension <- get_and_sanitize("PHENO_FILE_NAME_NO_EXTENSION")

args <- commandArgs(trailingOnly = TRUE)

PRSice_R <- gsub("\\\\", "", args[1])
PRSice_exec <- gsub("\\\\", "", args[2])
eigenvec_file <- gsub("\\\\", "", args[3])
snplist_file <- gsub("\\\\", "", args[4])
ibbc_analysis_prefix <- gsub("\\\\", "", args[5])

if (1 == 0) {

    setwd("~/IBBC/")
    R_dir <- "~/Dropbox/22Q11/R/"
    external_dir <- "external/"
    iBBC <- "external/22Qonly_INDM2_Sex_Mind001_SNPqc"
    ibbc_freq_file <- "2018_06_25/22Qonly_INDM2_Sex_Mind001_SNPqc.frqx"
    results_dir <- "2018_11_28"
    PRSice_R <- "bin//PRSice.R"
    PRSice_exec <- "bin//PRSice_mac"
    eigenvec_file <- "2018_06_25//pca.allsnps.eigenvec"
    snplist_file <- "2018_06_25//snplist.ibbc.txt"
    setwd("~/IBBC/")
    ibbc_analysis_prefix <- "2018_06_25/22Qonly_INDM2_Sex_Mind001_SNPqc.good_snps"
    pheno_file_name_no_extension <- "iBBC_AIMIIdata_14June2018"

}

source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "pheno_info.R"))

run_PRSice_all_phenotypes(target_file = ibbc_analysis_prefix, rebuild = FALSE)

## print_log_counts()

## load up phenotypes, add back on
pheno_csv_file <- file.path(external_dir, paste0(pheno_file_name_no_extension, ".csv"))
pheno <- read.csv(pheno_csv_file)
eigenvec <- read.table(eigenvec_file, header = TRUE)
pheno <- pheno[match_iid_to_affy_ids(eigenvec[, "IID"], pheno), ]
pheno <- cbind(pheno, eigenvec)
pheno <- add_all_prs_to_pheno() ## has default pvalue cutoff, but otherwise can import

output_file <- file.path(results_dir, paste0(pheno_file_name_no_extension, ".withPRS.csv"))
message(output_file)
write.csv(
    pheno,
    file = output_file,
    row.names = FALSE
)

## first p-value is 3.200558e-07
## R2 is 0.08747491 
case <- array(NA, nrow(pheno))
case[pheno$group2018 == "Case_SSD"] <- 1
case[pheno$group2018 == "Control"] <- 0
m <- rms::lrm(case ~  PRS_PGC_2014_SCZ, data = pheno)
print(m$stats["R2"])
s <- summary(glm(case ~ PRS_PGC_2014_SCZ, data = pheno, family = binomial))
coefficients(s)

## third p-value is from new PRS
case <- array(NA, nrow(pheno))
case[pheno$group2018 == "Case_SSD"] <- 1
case[pheno$group2018 == "Control"] <- 0
m <- rms::lrm(case ~  PRS_PGC_2018_SCZ, data = pheno)
print(m$stats["R2"])
s <- summary(glm(case ~ PRS_PGC_2018_SCZ, data = pheno, family = binomial))
coefficients(s)

## second p-value is 1.793111e-08 Adjusted R-squared: 0.04305204
s <- summary(lm(pheno[, "FSIQ_Z_First"] ~ PRS_Neale_XXXX_fluid, data = pheno))
print(s$r.squared)
coefficients(s)

