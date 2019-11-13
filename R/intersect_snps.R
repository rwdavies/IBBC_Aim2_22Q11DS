get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

R_dir <- get_and_sanitize("R_DIR")
results_dir <- get_and_sanitize("RESULTS_DIR")
external_dir <- get_and_sanitize("EXTERNAL_DIR")
iBBC <- get_and_sanitize("IBBC_PLINK")
pheno_file_name_no_extension <- get_and_sanitize("PHENO_FILE_NAME_NO_EXTENSION")
maf_filter <- as.numeric(get_and_sanitize("IBBC_MAF_FILTER"))
info_thresh <- as.numeric(get_and_sanitize("GWAS_MIN_INFO"))

args <- commandArgs(trailingOnly = TRUE)

ibbc_freq_file <- gsub("\\\\", "", args[1])
intersect_out_file <- gsub("\\\\", "", args[2])
ibbc_out_file <- gsub("\\\\", "", args[3])

if (1 == 0) {

    setwd("~/iBBC/")    
    R_dir <- "~/Dropbox/22Q11/R/"
    external_dir <- "external/"
    iBBC <- "external/22Qonly_INDM2_Sex_Mind001_SNPqc"
    ibbc_freq_file <- "2018_06_22/22Qonly_INDM2_Sex_Mind001_SNPqc.frqx"
    maf_filter <- 0.10
    info_thresh <- 0.90
    intersect_out_file <- tempfile()
    
}

file <- file.path(R_dir, "functions.R")
if (file.exists(file) == FALSE) {
    stop("bad R_dir")
}
source(file)
source(file.path(R_dir, "pheno_info.R"))

## now, load SNPs using rsid / other matching
iBBC_snps <- data.table::fread(paste0(iBBC, ".bim"), data.table = FALSE)
iBBC_freq <- data.table::fread(ibbc_freq_file, data.table = FALSE)
N <- sum(iBBC_freq[1, c("C(HOM A1)", "C(HET)", "C(HOM A2)")])
iBBC_snps$af <- (iBBC_freq[, "C(HET)"] + 2 * iBBC_freq[, "C(HOM A2)"]) / 2 / N

if (sum(iBBC_freq[, "SNP"] != iBBC_snps[, 2]) > 0) {
    stop("bad assumption")
}
## first, filter on MAF and region
in_mhc <-
    (iBBC_snps[, 1] == 6) &
    (26000000 <= iBBC_snps[, 4]) &
    (iBBC_snps[, 4] <= 34000000)
in_22q11 <-
    (iBBC_snps[, 1] == 22) &
    (18820303 <= iBBC_snps[, 4]) &
    (iBBC_snps[, 4] <= 21489474)
keep <-
    ((maf_filter < iBBC_snps$af) & ( iBBC_snps$af < (1 - maf_filter))) &
    (in_mhc == FALSE) &
    (in_22q11 == FALSE)
## know how many lost
message(sum(keep), " / ", length(keep))
message(sum(! keep), " / ", length(keep))
iBBC_snps2 <- iBBC_snps[keep, ]

## write out this SNPlist
write.table(
    matrix(iBBC_snps2[, 2], ncol = 1),
    file = ibbc_out_file,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE
)


## do this only on first two studies, for lack of change
## length(pheno_info)
out <- lapply(1:2, function(i_pheno) {
    file <- pheno_info[[i_pheno]][["file"]]
    if (substr(file, nchar(file) - 2, nchar(file)) == ".gz") {
        file <- paste0("gunzip -c ", file)
    }
    data <- data.table::fread(file, data.table = FALSE)
    info <- pheno_info[[i_pheno]][["info"]]
    if (is.na(info) == FALSE) {
        print("remove SNPs with low info")
        print(pheno_info[[i_pheno]]$pretty_name)        
        keep <- data[, "info"] > info_thresh
        print(paste0("keep ", sum(keep), " / ", length(keep)))
        print(paste0("remove ", sum(!keep), " / ", length(keep)))
        data <- data[keep, ]
    }
    return(data)
})

## now, intersect iBBC with all subsequent GWAS
## chr:bp:(allele 1
int <- intersect(iBBC_snps2[, 2], intersect(out[[1]][, 2], out[[2]][, 2]))

## check alleles
w1 <- match(int, iBBC_snps2[, 2])
a1 <- iBBC_snps2[w1, 5]
b1 <- iBBC_snps2[w1, 6]

## keep those 
i_pheno <- 1
w2 <- match(int, out[[i_pheno]][, pheno_info[[i_pheno]]$snp])
a2 <- out[[i_pheno]][w2, pheno_info[[i_pheno]]$effect_allele]
b2 <- out[[i_pheno]][w2, pheno_info[[i_pheno]]$non_effect_allele]

i_pheno <- 2
w3 <- match(int, out[[i_pheno]][, pheno_info[[i_pheno]]$snp])
a3 <- out[[i_pheno]][w3, pheno_info[[i_pheno]]$effect_allele]
b3 <- out[[i_pheno]][w3, pheno_info[[i_pheno]]$non_effect_allele]

## keep those that do not have disagreable alleles
keep1 <-
    ((a1 == a2) & (b1 == b2)) |
    ((a1 == b2) & (b1 == a2)) 
keep2 <-
    ((a1 == a3) & (b1 == b3)) |
    ((a1 == b3) & (b1 == a3))

## write out snplist
write.table(
    matrix(int[keep1 & keep2], ncol = 1),
    file = intersect_out_file,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t",
    quote = FALSE
)

check_af_between_iBBC_and_neale()

quit()






## remove ambiguous variants (A/T, C/G)?
table(a1, b1)
table(a2, b2) ## plenty?
table(a3, b3)

## accept those that
table(paste0(a1, "-", b1), paste0(a3, "-", b3))

iBBC_snps2[w1[1:10], ]
out[[2]][w3, ][1:10, ]

b1 <- iBBC_snps[match(int, iBBC_snps[, 2]), 6]


iBBC_snps2
