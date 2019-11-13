get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

R_dir <- get_and_sanitize("R_DIR")
results_dir <- get_and_sanitize("RESULTS_DIR")
external_dir <- get_and_sanitize("EXTERNAL_DIR")
iBBC <- get_and_sanitize("IBBC")
pheno_file_name_no_extension <- get_and_sanitize("PHENO_FILE_NAME_NO_EXTENSION")

args <- commandArgs(trailingOnly = TRUE)

snplist <- gsub("\\\\", "", args[1])
eigenvec_file <- gsub("\\\\", "", args[2])

if (1 == 0) {

    setwd("~/iBBC/")    
    R_dir <- "~/Dropbox/22Q11/R/"
    external_dir <- "external/"
    results_dir <- "~/iBBC/2018_06_22/"
    iBBC <- "external/22Qonly_INDM2_Sex_Mind001_SNPqc"
    ibbc_freq_file <- "2018_06_22/22Qonly_INDM2_Sex_Mind001_SNPqc.frqx"
    snplist <- "~/iBBC/2018_06_22/snplist.intersected.LDpruned.prune.in"
    pheno_file_name_no_extension <- "iBBC_AIMIIdata_14June2018"
    eigenvec_file <- "~/iBBC/2018_06_22/pca1.eigenvec"    
    
}

source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "pheno_info.R"))

pheno_csv_file <- file.path(external_dir, paste0(pheno_file_name_no_extension, ".csv"))
pheno <- read.csv(pheno_csv_file)
eigenvec <- read.table(eigenvec_file, header = TRUE)
pheno <- pheno[match_iid_to_affy_ids(eigenvec[, "IID"], pheno), ]
pheno <- cbind(pheno, eigenvec)

snps_to_keep <- as.character(read.table(snplist)[, 1])

genos <- BEDMatrix::BEDMatrix("~/iBBC/external/22Qonly_INDM2_Sex_Mind001_SNPqc.bed")
bim <- data.table::fread("~/iBBC/external/22Qonly_INDM2_Sex_Mind001_SNPqc.bim", data.table = FALSE)
## subset down
keep <- match(snps_to_keep, substr(colnames(genos), 1, nchar(colnames(genos)) - 2))
genos <- genos[, keep]
bim <- bim[keep, ]
## fill in missing with average
##af <- colSums(genos, na.rm = TRUE) / 2 / colSums(is.na(geno) == FALSE)
##for(i_snp in 1:ncol(genos)) {
## }

out <- lapply(1:length(pheno_info), function(i_pheno) {
    file <- pheno_info[[i_pheno]][["file"]]
    if (substr(file, nchar(file) - 2, nchar(file)) == ".gz") {
        file <- paste0("gunzip -c ", file)
    }
    data <- data.table::fread(file, data.table = FALSE)
    data <- data[match(snps_to_keep, data[, pheno_info[[i_pheno]]$snp]), ]
    return(data)
})


remove_ambiguous <- FALSE
for (i_pheno in 1:length(pheno_info)) {

    data <- out[[i_pheno]]
    ## get beta, match alleles
    if (pheno_info[[i_pheno]]$base_phenotype_is_binary) {
        beta <- log(data[, pheno_info[[i_pheno]]$or])
    } else {
        beta <- data[, pheno_info[[i_pheno]]$beta]
    }
    to_flip <- data[, pheno_info[[i_pheno]]$non_effect_allele] != bim[, 6]
    not_to_flip <- data[, pheno_info[[i_pheno]]$non_effect_allele] == bim[, 6]    
    #print(table(paste0(
    #    data[, pheno_info[[i_pheno]]$non_effect_allele], "-",
    #    data[, pheno_info[[i_pheno]]$effect_allele]),
    #    paste0(bim[, 5], "-", bim[, 6])))

    print(mean(beta))
    print(mean(beta[to_flip]))
    beta[to_flip] <- -beta[to_flip]
    print(mean(beta[to_flip]))    

    ## keep only those that meet p-value threshold
    print(pretty_name)
    pval_cutoff <- pheno_info[[i_pheno]]$pval_cutoff
    keep <- data[, pheno_info[[i_pheno]]$pvalue] < pval_cutoff        
    print(paste0(pretty_name, "- keep ", sum(keep), " / ", length(keep), " p-values meeting threshold"))

    if (remove_ambiguous) {
        remove <-
            (bim[, 5] == "A" & bim[, 6] == "T") |
            (bim[, 5] == "T" & bim[, 6] == "A") |
            (bim[, 5] == "C" & bim[, 6] == "G") |
            (bim[, 5] == "G" & bim[, 6] == "C")
        keep[remove] <- FALSE
        print(paste0(pretty_name, "-keep ", sum(keep), " / ", length(keep), " p-values meeting threshold"))        
    }

    prs <- genos[, keep] %*% matrix(beta[keep], ncol = 1)
    ## standardize
    prs <- (prs - mean(prs)) / sd(prs)
    pheno[, paste0("PRS_", names(pheno_info)[i_pheno])] <- prs[match(rownames(genos)    , paste0(pheno[, "FID"], "_", pheno[, "IID"]))]
    
}

output_file <- file.path(results_dir, paste0(pheno_file_name_no_extension, ".withPRS.csv"))
message(output_file)
write.csv(
    pheno,
    file = output_file,
    row.names = FALSE
)

## tests / checks

