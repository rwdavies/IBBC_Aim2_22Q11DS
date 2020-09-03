get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

R_dir <- get_and_sanitize("R_DIR")
results_dir <- get_and_sanitize("RESULTS_DIR")
pheno_file_name_no_extension <- get_and_sanitize("PHENO_FILE_NAME_NO_EXTENSION")

args <- commandArgs(trailingOnly = TRUE)

pheno_file_with_prs <- gsub("\\\\", "", args[1])
aim2A_hist_filename <- gsub("\\\\", "", args[2])
aim2A_main_filename <- gsub("\\\\", "", args[3])
aim2A_r2_filename <- gsub("\\\\", "", args[4])
aim2B_hist_filename <- gsub("\\\\", "", args[5])
aim2B_main_filename <- gsub("\\\\", "", args[6])
aim2X_tableSCZ_filename <- gsub("\\\\", "", args[7])
aim2X_tableIQ_filename <- gsub("\\\\", "", args[8])
clozuk_overlap <- gsub("\\\\", "", args[9])

if (isTRUE(as.logical(Sys.getenv("MANUAL_ANALYSIS_SWITCH")))) {

    R_dir <- "~/proj/IBBC_Aim2_22Q11DS/R/"
    results_dir <- file.path("~/IBBC/", "2018_11_28")
    external_dir <- file.path("~/IBBC/", "external")
    pheno_file_name_no_extension <- "iBBC_AIMIIdata_14June2018"
    pheno_file_with_prs <- file.path(results_dir, paste0(pheno_file_name_no_extension, ".withPRS.csv"))
    aim2A_hist_filename <- file.path(results_dir, "aim2A.hist.pdf")    
    aim2A_main_filename <- file.path(results_dir, "aim2A.main.pdf")
    aim2A_r2_filename <- file.path(results_dir, "aim2A.r2.pdf")
    aim2B_hist_filename <- file.path(results_dir, "aim2B.hist.pdf")        
    aim2B_main_scz_filename_pdf <- file.path(results_dir, "aim2B.main.scz.pdf")
    aim2B_main_scz_filename_png <- file.path(results_dir, "aim2B.main.scz.png")
    aim2B_main_iq_filename_pdf <- file.path(results_dir, "aim2B.main.iq.pdf")
    aim2B_main_iq_filename_png <- file.path(results_dir, "aim2B.main.iq.png")
    aim2B_main_new_filename_pdf <- file.path(results_dir, "aim2B.main.new.pdf")
    aim2B_main_new_filename_png <- file.path(results_dir, "aim2B.main.new.png")
    aim2B_main_old_filename_pdf <- file.path(results_dir, "aim2B.main.old.pdf")
    aim2B_main_old_filename_png <- file.path(results_dir, "aim2B.main.old.png")
    aim2X_tableSCZ_filename <- file.path(results_dir, "aim2X.tableSCZ.csv")
    aim2X_tableIQ_filename <- file.path(results_dir, "aim2X.tableIQ.csv")
    clozuk_overlap <- file.path(external_dir, "List_samples_Overlapping_with_CLOZUK.txt")

}

## library("DescTools")
library("VGAM")
source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "analysis_functions.R"))

pheno <- read.csv(pheno_file_with_prs)
## remove CLOZUK overlap
overlap_samples <- as.character(read.table(clozuk_overlap)[, 1])
pheno <- pheno[-match(overlap_samples, pheno[, "IID"]), ]
pheno$binary_VIQ_decline <- pheno[, "VIQ_deltazFL"] < (-1.0)
pheno$binary_FSIQ_decline <- pheno[, "FSIQ_deltazFL"] < (-1.0)
##
pheno$binary_SCZ_vs_merged <- array(NA, nrow(pheno))
pheno$binary_SCZ_vs_merged[(pheno[, "group2018"] == "Case_SSD") ] <- TRUE
pheno$binary_SCZ_vs_merged[(pheno[, "group2018"] == "Control") | (pheno[, "group2018"] == "PutativeControl")] <- FALSE
## 
pheno$binary_SCZ_vs_definite <- array(NA, nrow(pheno))
pheno$binary_SCZ_vs_definite[(pheno[, "group2018"] == "Case_SSD") ] <- TRUE
pheno$binary_SCZ_vs_definite[pheno[, "group2018"] == "Control"] <- FALSE
##
pheno$binary_sub_vs_merged <- array(NA, nrow(pheno))
pheno$binary_sub_vs_merged[(pheno[, "group2018"] == "PutativeSubthreshold") ] <- TRUE
pheno$binary_sub_vs_merged[(pheno[, "group2018"] == "Control") | (pheno[, "group2018"] == "PutativeControl")] <- FALSE
## make extra as well
pheno$binary_sub_vs_definite <- array(NA, nrow(pheno))
pheno$binary_sub_vs_definite[(pheno[, "group2018"] == "PutativeSubthreshold") ] <- TRUE
pheno$binary_sub_vs_definite[pheno[, "group2018"] == "Control"] <- FALSE
pheno$binary_sub_vs_definite[pheno[, "group2018"] == "PutativeControl"] <- NA

## write.table(pheno[, c("IID", "X", "genomics_id", "site_id", "group2018")], file = "~/IBBC/2018_11_28/992.samples.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)

get_sample_characteristics(
    output_csv = file.path(results_dir, "sample_characteristics.csv"),
    pheno = pheno
)
##pheno$binary_VIQ_decline <- as.integer(pheno[, "VIQ_deltazFL"] < -1.0)
##pheno$binary_FSIQ_decline <- as.integer(pheno[, "FSIQ_deltazFL"] < -1.0)

## I think this is just for certain histograms
groups_to_plot <- c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control", "Unknown", "Case_AffectivePsychosis")
phenoS2 <- make_pheno_with_ordered_group2018(pheno)

aim2A_results <- get_aim2A_results()     


source(file.path(R_dir, "analysis_functions.R"))
aim2A_plot_hist(aim2A_hist_filename)
for(type in c("pdf", "eps")) {
    aim2A_plot_groups(aim2A_main_filename, type = type)
}
aim2A_plot_r2(aim2A_r2_filename)

## check average levels of PRS for main text
sapply(c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control"), function(who) {
    mean(phenoS2[phenoS2[, "group2018"] == who, "PRS_PGC_2018_SCZ"])
})





##
## specifically for reviewer comment
## 
phenotypes <- c("binary_SCZ_vs_merged", "binary_sub_vs_definite", "FSIQ_Z_First", "binary_VIQ_decline")
names(phenotypes) <- c("SSD", "Subthreshold psychosis vs definite", "Baseline FSIQ", "Binary VIQ decline")
aim2B_results <- get_aim2B_results(phenotypes)
aim2B_plot_groups(paste0(aim2B_main_scz_filename_pdf, ".reviewer1.pdf"), "pdf", plot_type = "scz", add_beta = TRUE)
phenotypes <- c("binary_SCZ_vs_merged", "binary_sub_vs_merged", "FSIQ_Z_First", "binary_VIQ_decline")
names(phenotypes) <- c("SSD", "Subthreshold psychosis", "Baseline FSIQ", "Binary VIQ decline")
aim2B_results <- get_aim2B_results(phenotypes)
aim2B_plot_groups(paste0(aim2B_main_scz_filename_pdf, ".reviewer2.pdf"), "pdf", plot_type = "scz", add_beta = TRUE)

## so seems the same, likely from gaussian
coefficients(summary(glm(formula = binary_sub_vs_merged ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = pheno, family = gaussian)))["PRS_PGC_2014_SCZ", ]
coefficients(summary(glm(formula = binary_sub_vs_definite ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = pheno, family = gaussian)))["PRS_PGC_2014_SCZ", ]
coefficients(summary(glm(formula = binary_sub_vs_merged ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = pheno, family = binomial)))["PRS_PGC_2014_SCZ", ]
coefficients(summary(glm(formula = binary_sub_vs_definite ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = pheno, family = binomial)))["PRS_PGC_2014_SCZ", ]

##
## formulas <- c("binary_sub_vs_definite ~ PRS_PGC_2014_SCZ + sex + PC1 + PC2 + PC3 + PC4 + PC5", "binary_sub_vs_definite ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
## par(mfrow = c(1, 2))
## for(i in 1:2) {
##     s1 <- rms::lrm(formula = as.formula(formulas[i]), data = pheno)
##     y1 <- stats::predict(object = s1, newdata = pheno[w == TRUE & !is.na(w), ])
##     y2 <- stats::predict(object = s1, newdata = pheno[w == FALSE & !is.na(w), ])
##     plot(c(y1, y2), col = c(rep("blue", length(y1)), rep("red", length(y2))), main = c("Without maxassessmentage", "With maxassessmentage")[i], ylab = "y_i")
##     if (i == 2) legend("topright", col = c("blue", "red"), c("binary subthreshold", "definite control"), lwd = 2)
## }

## check power for hypothetical example
if (1 == 0) {
    
    set.seed(101)
    f <- function(mult, n) {
        x <- sapply(1:1000, function(iRep) {
            x1 <- rnorm(n, 0, sd = sqrt(0.90))
            x2 <- rnorm(n, 0, sd = sqrt(0.05))
            p <- 1 / (1 + exp(-(mult * x1 + x2)))
            z <- as.integer(runif(n) < p)
            if (iRep == 1) {
                print(rms::lrm(formula = z ~ x1 + x2, data = data.frame(z, x1, x2))$stats["R2"])
                print(cor(mult * x1 + x2, x1) ** 2)                
                print(cor(mult * x1 + x2, x2) ** 2)
            }
            b1 <- glm(formula = z ~ x1 + x2, data = data.frame(z, x1, x2), family = binomial)
            b2 <- glm(formula = z ~ x2, data = data.frame(z, x1, x2), family = binomial)            
            return(c(coefficients(summary(b1))[3, 4], coefficients(summary(b2))[2, 4]))
        })
        print(table(x[1, ] < 0.05))
        print(table(x[2, ] < 0.05))        
    }
    f(1, 500) ## 55% power
    f(6, 500) ## 22% power
    
}



## for ania, n=720 for baseline fsiq bit
## formula1 <- as.formula(paste0("case ~ ", prs, " + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5"))



keep <- rowSums(is.na(pheno[, c("FSIQ_Z_First", "PRS_PGC_2014_SCZ", "maxassessmentage", "sex", "PC1", "PC2", "PC3", "PC4", "PC5")])) == 0
s1 <- summary(lm(data = pheno, formula =  as.formula("FSIQ_Z_First ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")))
s2 <- summary(lm(data = pheno[keep, ], formula =  as.formula("FSIQ_Z_First ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")))
table(keep) ## 720
write.table(pheno[keep, "IID"], file = file.path(results_dir, "ania_FSIQ_Z_First_720_2020_05_19B.csv"), row.names = FALSE, col.names =TRUE, sep = ",", quote = TRUE)
## visual inspection the same
keep <- rowSums(is.na(pheno[, c("binary_VIQ_decline", "PRS_PGC_2014_SCZ", "maxassessmentage", "sex", "PC1", "PC2", "PC3", "PC4", "PC5")])) == 0
s1 <- summary(lm(data = pheno, formula =  as.formula("binary_VIQ_decline ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")))
s2 <- summary(lm(data = pheno[keep, ], formula =  as.formula("binary_VIQ_decline ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")))
table(keep) ## 396
write.table(pheno[keep, "IID"], file = file.path(results_dir, "ania_binary_viq_decline_396_2020_05_19B.csv"), row.names = FALSE, col.names =TRUE, sep = ",", quote = TRUE)
##
write.table(pheno, file = file.path(results_dir, "ania_pheno_allB.csv"), row.names = FALSE, col.names =TRUE, sep = ",", quote = TRUE)




##
## r2 and p-value for definite controls only
##
phenotypes <- c("binary_SCZ_vs_merged", "binary_SCZ_vs_definite", "binary_sub_vs_merged", "binary_sub_vs_definite")
names(phenotypes) <- c("SSD vs merged", "SSD vs definite", "Sub vs merged", "Sub vs definite")
aim2B_results <- get_aim2B_results(phenotypes)
aim2B_plot_groups(paste0(aim2B_main_scz_filename_pdf, ".definite.pdf"), "pdf", plot_type = "scz", add_beta = TRUE)



##
## normal progression here
## 
phenotypes <- c("binary_SCZ_vs_merged", "binary_sub_vs_merged", "FSIQ_Z_First", "binary_VIQ_decline")
names(phenotypes) <- c("SSD", "Subthreshold psychosis", "Baseline FSIQ", "Binary VIQ decline")
aim2B_results <- get_aim2B_results(phenotypes)
aim2B_plot_hist(aim2B_hist_filename)

##
source(file.path(R_dir, "analysis_functions.R"))
mc <- TRUE
for(add_beta in c(FALSE, TRUE)) {
    ## what are these - what phenotypes
    aim2B_plot_groups(aim2B_main_scz_filename_pdf, "pdf", plot_type = "scz", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(aim2B_main_scz_filename_png, "png", plot_type = "scz", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(aim2B_main_iq_filename_pdf, "pdf", plot_type = "iq", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(aim2B_main_iq_filename_png, "png", plot_type = "iq", add_beta = add_beta, mc = mc)
    ## now do subsets - very similar, just slightly different names I think
    phenotypes <- c("binary_SCZ_vs_merged", "FSIQ_Z_First")
    names(phenotypes) <- c("SSD", "Baseline FSIQ")
    aim2B_plot_groups(aim2B_main_old_filename_pdf, "pdf", plot_type = "both", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(aim2B_main_old_filename_png, "png", plot_type = "both", add_beta = add_beta, mc = mc)
    ## add for old? try it out?
    aim2B_plot_groups(filename = gsub("old", "old.special", aim2B_main_old_filename_pdf), what = "pdf", plot_type = "special_old", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(filename = gsub("old", "old.special", aim2B_main_old_filename_png), what = "png", plot_type = "special_old", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(filename = gsub("old", "old.special", gsub("png", "eps", aim2B_main_old_filename_png)), what = "eps", plot_type = "special_old", add_beta = add_beta, mc = mc)
    phenotypes <- c("binary_sub_vs_merged", "binary_VIQ_decline")
    names(phenotypes) <- c("Subthreshold psychosis", "VIQ decline")
    aim2B_plot_groups(aim2B_main_new_filename_pdf, "pdf", plot_type = "both", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(aim2B_main_new_filename_png, "png", plot_type = "both", add_beta = add_beta, mc = mc)
    ## do both but do in combined way
    aim2B_plot_groups(filename = gsub("new", "new.special", aim2B_main_new_filename_pdf), what = "pdf", plot_type = "special", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(filename = gsub("new", "new.special", aim2B_main_new_filename_png), what = "png", plot_type = "special", add_beta = add_beta, mc = mc)
    aim2B_plot_groups(filename = gsub("new", "new.special", gsub("png", "eps", aim2B_main_new_filename_png)), what = "eps", plot_type = "special", add_beta = add_beta, mc = mc)
}




add_beta <- FALSE; mc <- TRUE;    aim2B_plot_groups(aim2B_main_new_filename_pdf, "pdf", plot_type = "both", add_beta = add_beta, mc = mc)


## alternatively, do quads of two sets of phenotypes - known and unknown


1

source(file.path(R_dir, "analysis_functions.R"))
look_at_scz_quantile(aim2X_tableSCZ_filename, output_plot = TRUE)
look_at_iq_quantile(aim2X_tableIQ_filename, output_plot = TRUE)

for(type in c("pdf", "eps")) {
    filename <- file.path(results_dir, "cutoffPSile_both.pdf")
    if (type == "pdf") {
        pdf(filename, height = 4, width = 8)
    } else {
        filename <- gsub(".pdf", ".eps", filename)        
        cairo_ps(file = filename, onefile = FALSE, fallback_resolution = 600, height = 4, width = 8)
    }
    par(mfrow = c(1, 2))
    look_at_scz_quantile(aim2X_tableSCZ_filename, output_plot = FALSE)
    look_at_iq_quantile(aim2X_tableIQ_filename, output_plot = FALSE)
    dev.off()
}


plot_percentiles_vs_prevalence()  ## plot

## do mediation analysis
source(file.path(R_dir, "mediation.R"))
mediation()

## do supplemental figure 2, fraction of controls vs ps_sz
for_jacob_plot_fraction_of_controls_and_ps_scz()


## check 0.5 value, original choice
pheno$binary_VIQ_decline <- pheno[, "VIQ_deltazFL"] < (-0.5)
phenotypes <- c("binary_SCZ_vs_merged", "binary_sub_vs_merged", "FSIQ_Z_First", "binary_VIQ_decline")
names(phenotypes) <- c("SSD", "Subthreshold psychosis", "Baseline FSIQ", "Binary VIQ decline")
aim2B_results <- get_aim2B_results(phenotypes)
aim2B_plot_groups(paste0(aim2B_main_iq_filename_pdf, ".fornumbers.pdf"), "pdf", plot_type = "scz", add_beta = TRUE)




## quantitative phenotype analysis



## mediation analysis



1




ipw_Y <- function(x) {
    D <- 0.78 ## effecti
    B <- 17
    C <- 2 * ((0.80 + 0.95) / 2 - D)
    A <- log(C / (0.80 - D) - 1) / ((25 - 17))
    D + C / (1 + exp(-A * (x - B)))
}

sigmoid <- function() {
    ## fit sigmoid?
    ##
    ipw_Y(17) == (0.8 + 0.95) / 2
    ipw_Y(9) == (0.8)
    ipw_Y(25) == (0.95)
    ##
    x_range <- seq(1, 50, 1)
    y <- sapply(x_range, function(x) {
        ipw_Y(x)
    })
    ## assumed probability of being a control
    plot(x_range, y, type = "l", main = "P(control)")
    abline(v = 17, col = "red")
    abline(h = 0.875, col = "purple")    
    abline(v = 9, col = "red")
    abline(h = 0.80, col = "purple")    
    abline(v = 25, col = "red")
    abline(h = 0.95, col = "purple")        
}

test_ipw <- function(pheno, age_var = "maxassessmentage", prs_var = "PRS_PGC_2014_SCZ") {
    ## among putative controls, age vs PRS
    who <- pheno[, "group2018"] == "PutativeControl"
    ##plot(pheno[who, age_var], pheno[who, prs_var])
    age <- pheno[who, age_var]
    ##
    p_cont <- ipw_Y(age)
    p_case <- 1 - ipw_Y(age)
    prs <- pheno[who, prs_var]
    ## linear trend
    s1 <- summary(lm(prs ~ age))
    s2 <- summary(lm(prs ~ p_case))
    print(s1)
    print(s2)
    plot_trend <- function(s, col, z) {
        i <- coefficients(s)[1, 1]
        b <- coefficients(s)[2, 1]
        lines(age, i + b * z, col = col)
    }
    ## make plots?
    plot(age, prs, col = rgb(red = 0, green = 0, blue = 0, alpha = min(c(1, 1000 / (nrow(pheno))))), pch = 16)
    plot_trend(s1, col = "red", z = age)
    plot_trend(s2, col = "blue", z = 1 - Y(age))
    ## add linear trendline
}

## png("~/temp.png", height = 5, width = 15, res = 300, units = "in")
## par(mfrow = c(1, 3))
## sigmoid()
## test_ipw(pheno, "age", "full_scz") ## N = 150K, simulated
## test_ipw(phenoS2, "maxassessmentage", "PRS_PGC_2014_SCZ")
## dev.off()


quit()

table(pheno[, "group2018"], pheno[, "maxassessmentage"] >= 25)
pheno[which(pheno[, "maxassessmentage"] >= 25 & pheno[, "group2018"] == "PutativeSubthreshold"), c("IID", "maxassessmentage", "site_id", "group2018")]
##pheno[which(pheno[, "maxassessmentage"] < 25 & pheno[, "group2018"] == "Control"), c("IID", "maxassessmentage", "group2018")]



par(mfrow = c(2, 2))
for(group in c("Case_SSD", "Control","PutativeControl","PutativeSubthreshold" )) {
    w <- which(phenoS2[, "group2018"] == group)
    plot(sort(phenoS2[w, "maxassessmentage"]), main = group, ylim = c(0, 80))
}

groups <- c("Case_SSD", "Control","PutativeControl","PutativeSubthreshold" )
cutoffs <- seq(0, 125, by = 12.5)
age <- phenoS2[, "maxassessmentage"]
out <- t(sapply(1:(length(cutoffs) - 1), function(i) {
    w <- cutoffs[i] <= age & age < cutoffs[i + 1]
    x <- as.character(phenoS2[w, "group2018"])
    c <- NULL
    for(g in groups) {
        c <- c(c, sum(x == g, na.rm = TRUE))
        names(c)[length(c)] <- g
    }
    c <- c / sum(c)
    return(c)
}))
rownames(out) <- paste0(
    cutoffs[-length(cutoffs)],
    "_",
    cutoffs[-1]
)

pdf("~/temp.pdf", height = 5, width = 10)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
barplot(t(out), col = cbPalette[1:4], legend = groups, width = 1)
dev.off()
