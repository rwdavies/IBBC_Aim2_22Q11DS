##
    ## from analyze_prs.R
    ##
    R_dir <- "~/proj/IBBC_Aim2_22Q11DS/R/"
    source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "analysis_functions.R"))
    source(file.path(R_dir, "simulate_functions.R"))
    results_dir <- file.path("~/IBBC/", "2018_11_28")
    pheno_file_name_no_extension <- "iBBC_AIMIIdata_14June2018"
    pheno_file_with_prs <- file.path(results_dir, paste0(pheno_file_name_no_extension, ".withPRS.csv"))
pheno <- read.csv(pheno_file_with_prs)

## remove clozuk overlap
external_dir <- file.path("~/IBBC/", "external")
clozuk_overlap <- file.path(external_dir, "List_samples_Overlapping_with_CLOZUK.txt")
overlap_samples <- as.character(read.table(clozuk_overlap)[, 1])
pheno <- pheno[-match(overlap_samples, pheno[, "IID"]), ]

    groups_to_plot <- c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control", "Unknown", "Case_AffectivePsychosis")
    phenoS2 <- make_pheno_with_ordered_group2018(pheno)
    
    ## add in phenotype, examine
    phenoS3 <- phenoS2
    phenoS3$case <- NA
    phenoS3$case[phenoS2[, "group2018"] == "PutativeSubthreshold"] <- 1
    phenoS3$case[phenoS2[, "group2018"] == "PutativeControl"] <- 0
    phenoS3$case[phenoS2[, "group2018"] == "Control"] <- 0    
    formula1 <- as.formula("case ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
    results <- summary(glm(data = phenoS3, formula1, family = binomial))
    coefficients(results)["PRS_PGC_2014_SCZ", ]
    print(paste0("N = ", sum(results$df[1:2])))

    ## merge in mood here
    qsp <- read.csv("~/IBBC/external/QSP_from_jacob_2019_01_15.csv")
    par(mfrow = c(1, 2))
    qsp <- qsp[qsp[, "SIPS_entries"] >= 1, ]
    hist(qsp[, "SIPS1"])
    ## tiny normal
    ## hist(qsp[, "transformed_SIPS1"], breaks = 20)
    both <- merge(qsp, phenoS3, by = "genomics_id")
    ## OK, so remember for simulations
    table(both[, "group2018.x"])


    ##
    ## analysis here from "AIMII call minutes Jan 22. NEXT AIMII call: Tue Jan 29, 11 AM EST"
    ## 
    both[, "transformed_SIPS1"] <- qnorm(pexp(both[, "SIPS1"] + 0.5, rate = 0.2238))
    both[, "transformed_SIPS1"] <- (both[, "transformed_SIPS1"] - mean(both[, "transformed_SIPS1"])) / sd(both[, "transformed_SIPS1"])
    ## check looks good
    both$bin_sub <- NA
    both$bin_sub[both[, "group2018.x"] == "PutativeSubthreshold"] <- 1
    both$bin_sub[both[, "group2018.x"] == "PutativeControl"] <- 0
    ## double check histograms
    par(mfrow = c(1, 3))
    hist(both[, "transformed_SIPS1"])
    hist(both[both$bin_sub == 1, "transformed_SIPS1"])
    hist(both[both$bin_sub == 0, "transformed_SIPS1"])

    message("==================== quant vs PRS + binary =================")
    results <- lm(
        data = both,
        as.formula("transformed_SIPS1 ~ PRS_PGC_2014_SCZ + bin_sub + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
    )
    summary(results)
    
    message("==================== quant vs binary =================")    
    ## WITHOUT
    results <- lm(
        data = both,
        as.formula("transformed_SIPS1 ~ bin_sub + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
    )
    summary(results)
    
    ## N = 347
    ## p = 0.76
    ## r2 = 0.0001
    message("==================== binary vs PRS =================")
    results <- lm(
        data = both,
        as.formula("bin_sub ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
    )
    summary(results)
    
    message("==================== quant vs PRS =================")
    results <- lm(
        data = both,
        as.formula("transformed_SIPS1 ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
    )
    summary(results)




    
    message("===============")
    results <- lm(
        data = both[both$bin_sub == 1, ],
        as.formula("transformed_SIPS1 ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
    )
    summary(results)
    message("===============")    
    results <- lm(
        data = both[both$bin_sub == 0, ],
        as.formula("transformed_SIPS1 ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
    )
    summary(results)

    
    ## OK, do some power analyses here (or at least, do them in other window...)

    ## so 1) use quantitative phenotype
    ## 2) use quantiled phenotype
## pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
    a <- "PutativeSubthreshold"
    b <- "PutativeControl"
    g1 <- pheno[, "group2018"]
    g2 <- both[, "group2018.x"] ## REAL
    ## fit exponential?
    y <- both[(g2 == a) | (g2 == b), "SIPS1"]
    m <- table(y)
    x <- as.numeric(names(m)) + 0.5
    y <- as.numeric(as.character(m))
    y <- y / sum(y)
    fo3 <- y ~ 1/(b + x^c) # omit .lin parameter; plinear will add it automatically
    fo4 <- y ~ a * exp(-a * x)
    fm3 <- nls(fo3, data = data.frame(y = y, x = x), start = list(b = 1, c = 1), alg = "plinear")
    fm4 <- nls(fo4, data = data.frame(y = y, x = x), start = list(a = 0.5), alg = "plinear") ## this one
    plot(x = x, y = y)
    yf2 <- function(x) {
        0.8521 * 0.2238 * exp(-0.2238 * x)
    }
    lines(sapply(x, yf2), col = "blue")
    lines(0.9530 * dexp(x = x, rate = 0.2238), col = "purple")
    lines(dexp(x = x, rate = 0.2238), col = "orange")
    ## OK, so this is the distribution, basically

    ## so these are the probabilities of the group
    ## convert to average quantile
    pexp(q = x, rate = 0.2238)
    ## convert this to a normal
    y <- qnorm(pexp(q = x, rate = 0.2238))

    ## NOW UNDO


    
    ## TRANSFORM QUESTION -> NORMAL
    qsp[, "transformed_SIPS1"] <- qnorm(pexp(qsp[, "SIPS1"] + 0.5, rate = 0.2238))
    qsp[, "transformed_SIPS1"] <- (qsp[, "transformed_SIPS1"] - mean(qsp[, "transformed_SIPS1"])) / sd(qsp[, "transformed_SIPS1"])

    ## check OK
    library("VGAM")
    source(file.path(R_dir, "functions.R"))
    source(file.path(R_dir, "analysis_functions.R"))
    source(file.path(R_dir, "simulate_functions.R"))
    load(file = file.path(results_dir, "model_params.RData"))
    ##    pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno


    
    ## TRANSFORM NORMAL -> QUESTION
    ## SO to go from questions -> normal, do qnorm(pexp(q = x, rate = 0.2238))
## to go from normal -> questions, so
if (1 == 0) {
    png("~/IBBC/2018_11_28/quant_sub.png", height = 9, width = 12, units = "in",res = 300)
    par(mfrow = c(3, 4))
for(i in 1:3) {
    g <- NULL
        if (i == 1) {
            w <- (g == a) | (g == b); main = paste0(a, " | ", b)
            w2 <- (g2 == a) | (g2 == b)
        }
        if (i == 2) { w <- (g == a); main = a; w2 <- (g2 == a)}
        if (i == 3) { w <- (g == b); main = b; w2 <- (g2 == b)}
        ## SIPS QUESTION DISTRIBUTION (what I have naturally)
        ## SIPS NORMAL DISTRIBUTION
        ## QUANT QUESTION DISTRIBUTION
        ## QUANT NORMAL DISTRIBUTION (what I have naturally)
        cex.main <- 1.25
        hist(both[w2, "SIPS1"], main = paste0("Original: Question dist\n", main), cex.main = cex.main, col = "lightblue", xlim = c(0, 40), breaks = 20)
        hist(both[w2, "transformed_SIPS1"], main = paste0("Original: Transformed to normal\n", main), cex.main = cex.main, col = "orange", xlim = c(-3, 3), breaks = 21)
        p <- pheno[w, "sub"]
        p2 <- round(qexp(pnorm(p), rate = 0.2238))
        hist(p2, main = paste0("Simulated: Transformed to question dist\n", main), cex.main = cex.main, col = "lightblue", xlim = c(0, 40), breaks = 20)
        hist(p[w], main = paste0("Simulated: Normal dist\n", main), cex.main = cex.main, col = "orange", xlim = c(-3, 3), breaks = 21) 
    }
    ## 
    dev.off()
}




    pdf("~/IBBC/2018_11_28/quant_sub.ori_only.pdf", height = 9, width = 6)
    par(mfrow = c(3, 2))
    for(i in 1:3) {
        g <- NULL
        main <- c("Both subthreshold psychosis and\n putative control", "Subthreshold psychosis", "Putative control")[i]
        if (i == 1) {
            w <- (g == a) | (g == b); 
            w2 <- (g2 == a) | (g2 == b)
        }
        if (i == 2) { w <- (g == a); w2 <- (g2 == a)}
        if (i == 3) { w <- (g == b); w2 <- (g2 == b)}
        ## SIPS QUESTION DISTRIBUTION (what I have naturally)
        ## SIPS NORMAL DISTRIBUTION
        cex.main <- 1.25
        hist(both[w2, "SIPS1"], main = paste0("Original question dist\n", main), cex.main = cex.main, col = "lightblue", xlim = c(0, 40), breaks = 20, xlab = "SIPS1")
        hist(both[w2, "transformed_SIPS1"], main = paste0("Transformed to normal\n", main), cex.main = cex.main, col = "orange", xlim = c(-3, 3), breaks = 21, xlab = "transformed SIPS1")
    }
    ## 
    dev.off()


    ## jacob plot
    cbPalette1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    cbPalette2 <- sapply(cbPalette1, function(col) {
        x <- col2rgb(col, alpha = 0.5) / 255
        return(rgb(x[1], x[2], x[3], alpha = 0.25))
    })
    ## add slight jitter
    m <- mean(both[, "PRS_PGC_2014_SCZ"])
    col <- array(NA, nrow(both))
    col[both[, "group2018.x"] == "PutativeSubthreshold"] <- cbPalette2[2]
    col[both[, "group2018.x"] == "PutativeControl"] <- cbPalette2[3]
    pdf(file.path(results_dir, "sips1.ps_scz.pdf"), height = 10, width = 5)
    par(mfrow = c(2, 1))
    for(i_norm in 1:2) {
        if (i_norm == 1) {
            x <- both[, "SIPS1"]
            xlab <- "SIPS1"            
        } else {
            x <- qnorm(pexp(q = 0.5 + both[, "SIPS1"], rate = 0.2238))
            xlab <- "transformed SIPS1"
        }
        y <- both[, "PRS_PGC_2014_SCZ"]        
        add_jitter <- FALSE
        if (add_jitter) {
            s1 <- sd(x)
            x <- x + rnorm(n = nrow(both), mean = 0, sd = s1 / 10)        
            s2 <- sd(y)
            y <- y + rnorm(n = nrow(both), mean = 0, sd = s2 / 100)
        }
        plot(x, y, xlab = xlab, ylab = "PS_SZ", col = col, axes = FALSE, pch = 16, cex = 2)
        axis(1)
        axis(2)
        if (i_norm == 1) {
            legend("topright", c("Subthreshold psychosis", "Putative control"), col = cbPalette1[2:3], lwd = 4)
        }
        ## add trendline
        c <- coefficients(summary(lm(y ~ x)))
        abline(c[1, 1], c[2, 1])
    }
    dev.off()

