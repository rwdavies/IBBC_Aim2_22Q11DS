## manual re-formatting of Neale et al data

## see here for README
## https://docs.google.com/spreadsheets/d/1b3oGI2lUt57BcuHttWaZotQcI0-mBRPyZihz87Ms_No/edit#gid=275725118
## NOTE: Variants are listed in the form CHR:POS:REF:ALT, where the ALT allele is the effect allele in the model

## Update Nov 2018
## based on Aug 2018 analysis
## wget https://www.dropbox.com/s/daq8uxjwm56ihh3/20016_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 20016_raw.gwas.imputed_v3.both_sexes.tsv.bgz

get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

external_dir <- get_and_sanitize("EXTERNAL_DIR")

setwd(external_dir)
system("curl -L -O https://www.dropbox.com/s/shsiq0brkax886j/20016.assoc.tsv.gz?dl=0")
system("curl -L -O  https://www.dropbox.com/s/ehnp53rfqmp6xjg/variants.tsv?dl=0")
system("mv variants.tsv?dl=0 variants.tsv")
system("gzip -1 -f variants.tsv")

## re-format file 
file <- "20016.assoc.tsv.gz?dl=0"
data <- data.table::fread(paste0("gunzip -c ", file), data.table = FALSE)

file <- "variants.tsv.gz"
variants <- data.table::fread(paste0("gunzip -c ", file), data.table = FALSE)


## split up variant column - slow as molasses
first <- t(sapply(strsplit(as.character(data[, "variant"]), ":"), I))
data$chr <- first[, 1]
data$pos <- first[, 2]
data$non_effect_allele <- first[, 3]
data$effect_allele <- first[, 4]
data$info <- variants[match(variants[, "variant"], data[, "variant"]), "info"]
data$af <- variants[match(variants[, "variant"], data[, "variant"]), "AF"]


## dump back to disk
out_file <- file.path("fluid_intelligence.neale.tsv")
data.table::fwrite(
    data,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
)
system(paste0("gzip -1 -f ", out_file))

quit()


bim <- data.table::fread("/Users/robert davies/IBBC/2018_11_15/22Qonly_INDM2_Sex_Mind001_SNPqc.good_snps.bim", data.table = FALSE)
rownames(bim) <- paste(bim[, 1], bim[, 4], bim[, 5], bim[, 6], sep = ":")
both <- intersect(data[, "variant"], rownames(bim))
w1 <- match(both, rownames(bim))
w2 <- match(both, data[, "variant"])
table(
    paste0(bim[w1, 5], bim[w1, 6]),
    paste0(data[w2, "non_effect_allele"], data[w2, "effect_allele"])
)
           
