get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

R_dir <- get_and_sanitize("R_DIR")
external_dir <- get_and_sanitize("EXTERNAL_DIR")
pheno_file_name_no_extension <- get_and_sanitize("PHENO_FILE_NAME_NO_EXTENSION")
args <- commandArgs(trailingOnly = TRUE)

eigenvec_file <- gsub("\\\\", "", args[1])
pca_outliers <- gsub("\\\\", "", args[2])

if (isTRUE(as.logical(Sys.getenv("MANUAL_ANALYSIS_SWITCH")))) {

    R_dir <- "~/proj/IBBC_Aim2_22Q11DS/R/"
    external_dir <- "~/IBBC/external/"
    pheno_file_name_no_extension <- "iBBC_AIMIIdata_14June2018"
    eigenvec_file <- "~/IBBC/2018_11_28/pca.allsnps.eigenvec"
    pca_outliers <- "~/IBBC/2018_11_28/pca.allsnps.outliers"
    clozuk_overlap <- file.path(external_dir, "List_samples_Overlapping_with_CLOZUK.txt")    

}

eigenval_file <- gsub("eigenvec", "eigenval", eigenvec_file)

source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "pheno_info.R"))

cbPaletteR <- rep(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"), 100)
pchPaletteR <- rep(1:100, each = 8)

eigen <- read.table(eigenvec_file, header = TRUE)

## any difference?
pheno_csv_file <- file.path(external_dir, paste0(pheno_file_name_no_extension, ".csv"))
pheno <- read.csv(pheno_csv_file)
## 
eigenvec <- read.table(eigenvec_file, header = TRUE)
eigenval <- read.table(eigenval_file)
pheno <- pheno[match_iid_to_affy_ids(eigenvec[, "IID"], pheno), ]
pheno <- cbind(pheno, eigenvec)
pheno$site_id_colour <- cbPaletteR[match(pheno[, "site_id"], unique(pheno[, "site_id"]))]
pheno$site_id_pch <- pchPaletteR[match(pheno[, "site_id"], unique(pheno[, "site_id"]))]
##
## do removal now
##
overlap_samples <- as.character(read.table(clozuk_overlap)[, 1])
to_remove <- sapply(overlap_samples, function(x) grep(x, pheno[, "affy_id_1"]))
pheno <- pheno[-to_remove, ]



## visualize
w <- pheno[, "group2018"] %in% c("Case_SSD", "Control", "PutativeControl", "PutativeSubthreshold"); table(pheno[w, "site_id"], as.character(pheno[w, "group2018"]))

## eigenvalues
paste0(paste0(1:20, "=", round(eigenval[, 1], 3)), collapse = ', ')



pdf("~/IBBC/2018_11_28/pc1.fancy.pdf", height = 8, width = 6)
## for PC1, make single plot, use my own colours
par(mar = c(5, 0, 0, 0))
##par(oma = c(0, 3, 0, 0))
main <- NULL
## remove those with sites "DUPS"
pc_col <- "PC1"
phenoX <- pheno[order(pheno[, "site_id"], pheno[, pc_col]), ]
phenoX2 <- phenoX[!phenoX[, "site_id"] == "BBC_DUPS", ]
##
sites <- unique(phenoX2[, "site_id"])
##
y <- match(phenoX2[, "site_id"], sites) + rnorm(n = nrow(phenoX2), mean = 0, sd = 0.1)
col <- phenoX2$site_id_colour[match(sites, phenoX2[, "site_id"])]
names(col) <- sites
##y <- 1:nrow(phenoX2)
x <- phenoX2[, pc_col]
xlim <- c(min(x), max(x))
label_spot <- xlim[2] + 0.15 * diff(xlim)
xlim[2] <- xlim[2] + 0.25 * diff(xlim)
plot(x = x, y = y, col = phenoX2$site_id_colour, pch = phenoX2$site_id_pch, main = main, xlim = xlim, axes = FALSE, xlab = "PC1", ylab = "")
axis(1)
## labels
av_per_site <- sapply(sites, function(site) {
    mean(y[phenoX2[, "site_id"] == site])
})
names(av_per_site) <- sites
## add per-group average
for(site in sites) {
    a <- mean(x[(phenoX2[, "site_id"] == site)])
    segments(x0 = a, x1 = a, y0 = av_per_site[site] - 0.2, y1 = av_per_site[site] + 0.2, col = "black", lwd = 2)
}
text(x = label_spot, y = av_per_site, labels = sites, col = col)
## add spacing? or specific height?
dev.off()

print("NO OUTLIERS FROM VISUAL INSPECTION")
print("IF ON VISUAL INSPECTION, OUTLIERS, RE-WRITE THIS")
write.table(
    "",
    file = pca_outliers
)


quit()

## below used to be main
## but not needed for main analysis now


for(suffix in c("pdf", "png")) {
    image_open(filename = paste0(pca_outliers, ".several"), height = 40 * 0.6, width = 10 * 0.6, suffix = suffix)
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
}

for(suffix in c("pdf", "png")) {
    for(i_pc in 1:10) {
        pc_col <- paste0("PC", i_pc)
        pc_eigen <- eigenval[i_pc, 1]
        image_open(paste0(pca_outliers, ".", pc_col), height = 10, width = 20, suffix = suffix, res = 100)
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
        main <- paste0(pc_col, ", eigenval=", pc_eigen)
        plot(x = x, y = y, col = phenoX$site_id_colour, pch = phenoX$site_id_pch, main = main, xlab = "")
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
}

