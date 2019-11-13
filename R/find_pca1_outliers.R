get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

R_dir <- get_and_sanitize("R_DIR")
external_dir <- get_and_sanitize("EXTERNAL_DIR")
pheno_file_name_no_extension <- get_and_sanitize("PHENO_FILE_NAME_NO_EXTENSION")
args <- commandArgs(trailingOnly = TRUE)

pca1_eigenvec <- gsub("\\\\", "", args[1])
pca1_outliers <- gsub("\\\\", "", args[2])

if (1 == 0) {

    R_dir <- "~/Dropbox/22Q11/R/"
    external_dir <- "~/iBBC/external/"
    pheno_file_name_no_extension <- "iBBC_AIMIIdata_14June2018"
    pca1_eigenvec <- "~/iBBC/2018_06_22/pca1.eigenvec"
    pca1_outliers <- "~/iBBC/2018_06_22/pca1.eigenvec.outliers"
    
}

source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "pheno_info.R"))

cbPaletteR <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 100)
pchPaletteR <- rep(1:100, each = 8)

eigen <- read.table(pca1_eigenvec, header = TRUE)

## any difference?
pheno_csv_file <- file.path(external_dir, paste0(pheno_file_name_no_extension, ".csv"))
pheno <- read.csv(pheno_csv_file)
eigenvec <- read.table(eigenvec_file, header = TRUE)
eigenval <- read.table(eigenval_file)
pheno <- pheno[match_iid_to_affy_ids(eigenvec[, "IID"], pheno), ]
pheno <- cbind(pheno, eigenvec)
pheno$site_id_colour <- cbPaletteR[match(pheno[, "site_id"], unique(pheno[, "site_id"]))]
pheno$site_id_pch <- pchPaletteR[match(pheno[, "site_id"], unique(pheno[, "site_id"]))]

## visualize
w <- pheno[, "group2018"] %in% c("Case_SSD", "Control", "PutativeControl", "PutativeSubthreshold"); table(pheno[w, "site_id"], as.character(pheno[w, "group2018"]))

## eigenvalues
paste0(paste0(1:20, "=", round(eigenval[, 1], 3)), collapse = ', ')

pdf("~/iBBC/2018_06_22/pca.1.several.pdf", height = 40 * 0.6, width = 10 * 0.6)
par(mfrow = c(8, 2))
for(i in 1:8) {
    a <- paste0("PC", i)
    b <- paste0("PC", i + 1)
    c <- paste0("PC", i + 2)
    plot(x = pheno[, a], y = pheno[, b], xlab = a, ylab = b, pch = pheno[, "site_id_pch"], col = pheno[, "site_id_colour"])
    text(x = pheno[, a], y = pheno[, b], labels = 1:nrow(pheno), cex = 0.25)
    plot(x = pheno[, a], y = pheno[, c], xlab = a, ylab = c, pch = pheno[, "site_id_pch"], col = pheno[, "site_id_colour"])
    text(x = pheno[, a], y = pheno[, c], labels = 1:nrow(pheno), cex = 0.25)    
}
dev.off()

for(i_pc in 1:10) {
    pc_col <- paste0("PC", i_pc)
    pc_eigen <- eigenval[i_pc, 1]
    pdf(paste0("~/iBBC/2018_06_22/pca.1.", pc_col, ".pdf"), height = 10, width = 30)
    par(mfrow = c(2, 1))
    ## boxplot
    plot(pheno[, "site_id"], pheno[, pc_col])
    ## similar, but coloured
    phenoX <- pheno[order(pheno[, "site_id"], pheno[, pc_col]), ]
    x <- 1:nrow(phenoX)
    y <- phenoX[, pc_col]
    ylim <- c(
        min(min(y), -0.1),
        max(max(y), 0.1)
    )
    plot(x = x, y = y, col = phenoX$site_id_colour, pch = phenoX$site_id_pch, main = paste0(pc_col, ", eigenval=", pc_eigen))
    a <- unique(pheno[, "site_id"])
    for(i_a in 1:length(a)) {
        colL <- phenoX[which.max(phenoX$site_id == a[i_a]), "site_id_colour"]
        pchL <- phenoX[which.max(phenoX$site_id == a[i_a]), "site_id_pch"]
        xx <- mean(x[phenoX$site_id == a[i_a]])
        yy <- 0
        text(
            x = xx,
            y = yy,
            label = a[i_a],
            col = colL,
            srt = 90
        )
        ##    points(x = xx, y = yy, col = colL, pch = pchL)
    }
    dev.off()
}


