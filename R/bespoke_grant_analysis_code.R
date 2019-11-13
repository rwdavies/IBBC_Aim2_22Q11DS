###########################
## copied and pasted from analyze_prs.R
###########################
R_dir <- "~/Dropbox/22Q11/R/"
results_dir <- file.path("~/IBBC/", "2018_07_05")
pheno_file_name_no_extension <- "iBBC_AIMIIdata_14June2018"
pheno_file_with_prs <- file.path(results_dir, paste0(pheno_file_name_no_extension, ".withPRS.csv"))
aim2A_hist_filename <- file.path(results_dir, "aim2A.hist.pdf")    
aim2A_main_filename <- file.path(results_dir, "aim2A.main.pdf")
aim2A_r2_filename <- file.path(results_dir, "aim2A.r2.pdf")
aim2B_hist_filename <- file.path(results_dir, "aim2B.hist.pdf")        
aim2B_main_filename <- file.path(results_dir, "aim2B.main.pdf")    
source(file.path(R_dir, "analysis_functions.R"))
pheno <- read.csv(pheno_file_with_prs)
pheno$binary_VIQ_decline <- as.integer(pheno[, "VIQ_deltazFL"] > 0.5)
pheno$binary_FSIQ_decline <- as.integer(pheno[, "FSIQ_deltazFL"] > 0.5)
###########################


## 1. For the narrative of the grant it would be helpful to have the explained variance (pseudo R2) for the PRS_IQ on baseline IQ, with the latter as a dichotomous outcome (<70 vs >70). 
s <- summary(lm(VIQ_First < 70 ~ maxassessmentage + maxassessmentage ** 2 + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = pheno))
summary(lm(s$residuals ~ pheno[as.integer(names(s$residuals)), "PRS_Neale_XXXX_fluid"]))





## 2. If possible, can you generate a polygenic risk figure (for IQ, dichotomous) for different p-value thresholds, such as the one below? It will only serve the purpose of presenting pilot data for this grant. 
more <- read.table("~/IBBC/2018_07_05/Neale_XXXX_fluid_iBBC/PRSice.all.score", header = TRUE)
## generate levels
both <- merge(pheno, more, by = "IID")
## do same thing as above
default_model <- summary(lm(VIQ_First < 70 ~ maxassessmentage + maxassessmentage ** 2 + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = both))
pval_thresh <- c("X0.050000", "X0.100000", "X0.150000", "X0.200000", "X1.000000")
pval_thresh_clean <- c("0.05", "0.10", "0.15", "0.20", "1")
results <- sapply(pval_thresh, function(pval) {
    x <- summary(lm(default_model$residuals ~ both[as.integer(names(s$residuals)), pval]))
    c(
        adj.r.squared = x$adj.r.squared,
        p = x$coefficients[2, 4]
    )
})

## now make plot
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col <- colorRampPalette(c(cbPalette[2], cbPalette[3]))(ncol(results))
x <- seq(0, 1, length.out = ncol(results) + 1)
xleft <- x[-length(x)]
xright <- x[-1]
xav <- (xleft + xright )/ 2
ybottom <- rep(0, ncol(results))
ytop <- results["adj.r.squared", ]
pdf("~/IBBC/grant.pdf", height = 5, width = 5)
plot(x = 0, y = 0, xlim = c(0, 1), ylim = c(0, 1.5 * max(results["adj.r.squared", ])), axes = FALSE, xlab = "p-value thresholds", ylab = "Adjusted R2", col = "white", main = "Adjusted for age, age^2, sex & 5 PCs")
axis(2)
mtext(text = pval_thresh_clean, side = 1, at = xav)
rect(xleft, ybottom, xright, ytop, col = col, border = NULL)
## text above
text(
    x = xav,
    y = results["adj.r.squared", ],
    labels = paste0(
        "p = ", 
        formatC(results["p", ], format = "e", digits = 2)
    ),
    srt = 90,
    pos = 3,
    offset = 3
)
dev.off()




## 3. I really would like to chat with you about the feasibility of individual risk prediction using the PRS. For that purpose it would be interesting to divide the PRS distribution into different bins and calculate the OR for low IQ (<70) of the highest PRS bin (e.g. top 10% PRS score). See recent paper by Khera et al (attached).  
## As you may recall, I mentioned previously that the positive predictive value (PPV) of a test is highly influenced by the a priori prevalence of the disease phenotype. I would therefore expect a much more interesting (and potentially clinically relevant) PPV derived from the PRS in the 22q11DS population compared with the a PRS with comparable properties (R2, OR) in the general population. What it means is that the PRS for low IQ (or for Schizophrenia) may have clinical value (bc substantial PPV) whereas the PRS for low IQ in the gen pop is interesting but not relevant for risk prediction. 

## calculate OR for different bins, Khera style
default_model <- summary(lm(VIQ_First < 70 ~ maxassessmentage + maxassessmentage ** 2 + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = pheno))
x <- pheno[as.integer(names(default_model$residuals)), "PRS_Neale_XXXX_fluid"]
thresholds <- quantile(x, probs = c(0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50))
message("odds ratio versus remainder of population")
t(sapply(thresholds, function(t) {
    s <- summary(lm(default_model$residuals ~ (x < t)))
    or <- exp(s$coefficients[2, 1])
    p <- s$coefficients[2, 4]
    return(c(or = or, p = p))
}))



