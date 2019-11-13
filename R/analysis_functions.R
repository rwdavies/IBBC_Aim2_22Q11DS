make_PSile_plot <- function(
    out,
    out2,
    filename,
    groups,
    cutoffs,
    cutoffs2,
    altPrev,
    case_col,
    control_groupname,
    flip_sign,
    main,
    output_plot = TRUE
) {
    ##
    local_cbPalette <- cbPalette[c(4, 7)]
    ## get them here!
    f <- function(w) {
        return(c(
            sapply(out$results, function(x) x[w, case_col]),
            sapply(out2$results, function(x) x[w, case_col])
        ))
    }
    ppvs <- rbind(
        f("ppv"), f("ppv_lower"), f("ppv_upper"), f("ppv_altPrev"),
        f("tp"), f("fp")
    )
    b <- groups %in% c(case_col, control_groupname)
    a <- sum(groups == case_col, na.rm = TRUE) / sum(b, na.rm = TRUE)
    ppvs <- cbind(ppvs[, 1:3], c(rep(a, 3), altPrev, NA, NA), ppvs[, 4:6])
    ppvs <- rbind(ppvs, ppvs[5, ] + ppvs[6, ])
    rownames(ppvs) <- c(
        "ppv", "ppv_lower", "ppv_upper", "ppv_altPrev",
        "tp", "fp", "condP"
    )
    ppvs <- ppvs[, ncol(ppvs):1]
    if (flip_sign) {
        ## more manual than I'd like!
        colnames(ppvs) <- c(
            paste0(100 * cutoffs2[length(cutoffs2):1], "<"), #, "th\nile")
            "All",
            paste0("<", 100 * cutoffs[length(cutoffs):1])
        )
    } else {
        colnames(ppvs) <- c(
            paste0("<", 100 * (1 - cutoffs)), #, "th\nile"),
            "All",
            paste0(100 * (1 - cutoffs2), "<") #, "th\nile")
        )
        ##ppvs[, 1:length(cutoffs)] <- ppvs[, length(cutoffs):1]
        ##ppvs[, length(cutoffs) + 1:length(cutoffs2)] <- ppvs[, length(cutoffs) + length(cutoffs2):1]
    }
    if (output_plot) {
        png(filename, height = 5, width = 5, units = "in", res = 300)
    }
    xlim <- c(1, ncol(ppvs))
    ylim <- c(0, 0.8)
    x <- 1:ncol(ppvs)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, xlab = "", ylab = "Positive predictive value", col = "white", axes = FALSE, main = main)
    axis(1, labels = FALSE)
    axis(1, at = x, labels = colnames(ppvs), padj = 1, srt = 90)
    ## mtext(text = colnames(ppvs), at = x, padj = 3, las = 1.5, side = 1)
    mtext("Subgroups by different percentile PS cutoff", side=1, line=3)
    axis(2)
    ##
    for(i22 in c(1, 2)) {
        if (i22 == 1) {
            ppv <- ppvs[1, ]
            lower <- ppvs[2, ]
            upper <- ppvs[3, ]
        } else {
            ## general population here
            ppv <- ppvs[4, ]
            ## OK - try it out!
            A <- 1.96 * sqrt(ppv * (1 - ppv) / ppvs["condP", ])
            lower <- ppv - A
            upper <- ppv + A
            lower["All"] <- ppv["All"]
            upper["All"] <- ppv["All"]
            lower[lower < 0] <- 0
            print("-------general-------")
            print(
                rbind(
                    colnames(ppvs),
                    ppv,
                    upper,
                    lower
                )
            )
        }
        points(x = x, y = ppv, pch = 16, col = local_cbPalette[i22], cex = 0.5, type = "o")
        ## add ci here
        arrows(x0 = x, x1 = x, y0 = lower, y1 = upper, length = 0, col = local_cbPalette[i22])
        ## top bit
        w <- 0.05
        which <- (ppv != lower)
        x0 <- (x - w)[which]
        x1 <- (x + w)[which]
        for(j in 1:2) {
            if (j == 1) { y0 <- y1 <- upper[which] }
            if (j == 2) { y0 <- y1 <- lower[which] }
            arrows(x0 = x0, x1 = x1, y0 = y0, y1 = y1, length = 0, col = local_cbPalette[i22])
        }
        ## average values here
        abline(h = ppv["All"], col = local_cbPalette[i22], lty = 2)
        ##
    }
    legend("topleft", c("22q11.2DS", "General"), col = local_cbPalette, lwd = 2)
    ## average orange dotted line
    if (output_plot) {
        dev.off()
    }
}






## must be run interactively
## PseudoR2 requires global environment ARGH
get_aim2A_results_pretty <- function() {

    ## wrap in a loop to just get printout at the end
    for(i in 1:1) {
        ##
        ##
        ##
        for(
            what in
            c(
                "Case_SSD_vs_Control",
                "Case_SSD_vs_PutativeControl&Control",
                "PutativeSubthreshold_vs_PutativeControl&Control"
            )
        ) {
            message("--------")
            message(what)
            local_results <- aim2A_results[[what]]
            N <- length(summary(local_results[["with_age"]])$deviance.resid)
            phenoS <- local_results[["phenoS"]]
            x1 <- local_results[["with_age"]]
            x2 <- local_results[["with_age_without_prs"]]
            r2 <- PseudoR2(x1, "Nagelkerke") - PseudoR2(x2, "Nagelkerke")
            print(coefficients(summary(x1))[2, ])
            message(paste0("N = ", N))
            message(paste0("r2 = ", r2))
        }
        ##
        ## try a different threshold for binary VIQ decline and PRS_SCZ
        ##
        message("----------")
        message("look at different threshold for binary VIQ decline")
        pheno$binary_VIQ_decline_0.5 <- as.integer(pheno[, "VIQ_deltazFL"] < (-0.5))
        phenotypes <- c("binary_VIQ_decline_0.5")
        names(phenotypes) <- c("Binary VIQ decline alternate")
        aim2B_results_0.5 <- get_aim2B_results(phenotypes)
        ##
        prsA <- "scz"
        prs1 <- paste0("prs_", prsA)
        no_prs <- paste0(prsA, "_no_PRS")
        i_what <- 1
        x <- aim2B_results_0.5[[i_what]][[prs1]]
        print(coefficients(x)[2, ])
        x_noPRS <- aim2B_results_0.5[[i_what]][[no_prs]]
        ##
        r2A <- x$r.squared - x_noPRS$r.squared
        r2B <-
            (1 - x$deviance / x$null.deviance) -
            (1 - x_noPRS$deviance / x_noPRS$null.deviance)
        r2 <- c(r2A, r2B)
        print(r2)
        ##
        ## get N, p, r2 for subthreshold psychosis and PS_IQ
        ##
        aim2A_results_iq <- get_aim2A_results(prs = "PRS_Neale_XXXX_fluid")
        message("==================================")
        for(what in c("Case_SSD_vs_Control", "PutativeSubthreshold_vs_PutativeControl&Control")) {
            message("---- PS_IQ ----")
            message(what)
            local_results <- aim2A_results_iq[[what]]
            N <- length(summary(local_results[["with_age"]])$deviance.resid)
            phenoS <- local_results[["phenoS"]]
            x1 <- local_results[["with_age"]]
            x2 <- local_results[["with_age_without_prs"]]
            r2 <- PseudoR2(x1, "Nagelkerke") - PseudoR2(x2, "Nagelkerke")
            print(coefficients(summary(x1))[2, ])
            message(paste0("N = ", N))
            message(paste0("r2 = ", r2))
        }
        ##
        ## get average SCZ for different groups
        ##
        message("======== average PS_SZ for different groups ===================")
        for(
            who in
            c(
                "Case_SSD", "PutativeSubthreshold",
                "PutativeControl", "Control"
            )
        ){
            x <- pheno[pheno[, "group2018"] == who, "PRS_PGC_2014_SCZ"]
                print(paste0(
                    "mean(PS_SZ) = ",
                    mean(x, na.rm = TRUE),
                    ", N = ", sum(is.na(x) == FALSE)
                ))
        }
    }


}


get_sample_characteristics <- function(output_csv, pheno) {

    sf <- function(a, d1 = 1, d2 = 0) {
        n <- sum(is.na(a) == FALSE)
        paste0(
            round(mean(a, na.rm = TRUE), d1), " (",
            round(sd(a, na.rm = TRUE), d1), ") [",
            n, "], {",
            round(min(a, na.rm = TRUE), d2), ",",
            round(max(a, na.rm = TRUE), d2), "}"
        )
    }
    sfB <- function(a, d1 = 1, d2 = 0) {
        n <- sum(is.na(a) == FALSE)
        paste0(
            round(100 * sum(a == TRUE, na.rm = TRUE) / n, d1), "% (",
            n, ")"
        )
    }
    mood <- read.csv("~/IBBC/external/AIMII_mood_2018_01_22.csv")
    pheno$psy_mood <- mood[match(pheno[, "IID"], mood[, "IID"]), "psy_mood"]
    pheno$psy_mood <- c(NA, FALSE, TRUE, NA)[match(pheno$psy_mood, c("", "0", "1", "unk"))]
    ## make table1 of summary characteristics
    output_phenotypes <- c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control", "any", "FSIQ_Z_First", "binary_VIQ_decline")
    rm(results)
    linker <- c(
        percent_male = "sex",
        assessment_age = "maxassessmentage",
        first_FSIQ_age = "Age_FSIQ_First",
        last_FSIQ_age = "Age_FSIQ_Last",
        first_VIQ_age = "Age_VIQ_First",
        last_VIQ_age = "Age_VIQ_Last",
        percent_mood = "psy_mood",
        mood2 = "psy_mood",
        FSIQ = "FSIQ_Z_First",
        VIQ = "binary_VIQ_decline",
        BVIQ = "binary_VIQ_decline"
    )
    results <- sapply(
        1:length(output_phenotypes),
        function(i_who) {
            who <- output_phenotypes[i_who]
            if (i_who <= 4) {
                w <- pheno[, "group2018"] == who
            } else if (i_who == 5) {
                w <- pheno[, "group2018"] %in% c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control")
            } else {
                w <- !is.na(pheno[, who])
            }
            ## sex 2 = female, 1 = male
            to_out <- c(
                N = sum(w),
                percent_male = round(100 * sum(pheno[w, "sex"] == 1) / sum(w)),
                assessment_age = sf(pheno[w, linker["assessment_age"]]),
                first_FSIQ_age = sf(pheno[w, "Age_FSIQ_First"]),
                last_FSIQ_age = sf(pheno[w, "Age_FSIQ_Last"]),
                first_VIQ_age = sf(pheno[w, "Age_VIQ_First"]),
                last_VIQ_age = sf(pheno[w, "Age_VIQ_Last"]),
                percent_mood = round(100 * sum(pheno[w, "psy_mood"] == 1) / sum(w)),
                mood2 = sfB(pheno[w, "psy_mood"] == 1),
                FSIQ = sf(pheno[w, "FSIQ_Z_First"], d1 = 2),
                VIQ = sf(pheno[w, "binary_VIQ_decline"]),
                BVIQ = sfB(pheno[w, "binary_VIQ_decline"])
            )
            ##
            ## add more columns!
        }
    )
    ##results <- t(results)
    colnames(results) <- output_phenotypes
    ## add anova, for each column, using first 4 valuies
    results2 <- cbind(results, NA)
    colnames(results2)[ncol(results2)] <- "p"
    ## add in anova
    w <- pheno[, "group2018"] %in% c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control")
    for(p1 in names(linker)) {
        p2 <- linker[p1]
        y <- pheno[w, p2]
        res.aov <- aov(y ~ pheno[w, "group2018"], data = pheno)
        m <- summary(res.aov)[[1]]
        results2[p1, "p"] <- m[1, "Pr(>F)"]
    }
    results2 <- results2[
        c("N", "percent_male", "assessment_age", "FSIQ",
          "BVIQ", "mood2"), ]
    write.table(
        results2,
        file = output_csv,
        sep = ",",
        col.names = TRUE,
        row.names = TRUE
    )

    if (1 == 0) {

        ## get values, new table, for rest
        results2 <- sapply(
            phenotypes,
            function(phenotype) {
                w <- !is.na(pheno[, phenotype])
                c(
                    N = sum(w),
                    mean_age = mean(pheno[w, "maxassessmentage"], na.rm = TRUE),
                    min_age = min(pheno[w, "maxassessmentage"], na.rm = TRUE),
                    max_age = max(pheno[w, "maxassessmentage"], na.rm = TRUE)
                )
            }
        )

    }

}


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

make_pheno_with_ordered_group2018 <- function(pheno) {
    groups_to_plot <- c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control", "Unknown", "Case_AffectivePsychosis")
    phenoS2 <- pheno
    phenoS2 <- phenoS2[order(match(phenoS2[, "group2018"], groups_to_plot)), ]
    phenoS2 <- phenoS2[phenoS2[, "group2018"] %in% groups_to_plot[1:4], ]
    phenoS2[, "group2018"] <- factor(
        as.character(phenoS2[, "group2018"]),
        groups_to_plot[1:4]
    )
    return(phenoS2)
}

get_aim2A_results <- function(prs = "PRS_PGC_2014_SCZ") {
    vars1 <- list(
        c("Case_SSD"),
        c("PutativeSubthreshold"),
        c("Case_SSD"),
        c("PutativeControl"),
        c("Case_SSD", "PutativeSubthreshold"),
        c("PutativeSubthreshold"),
        c("PutativeSubthreshold"),
        c("Case_SSD")
    )
    vars2 <- list(
        c("Control"),
        c("Control"),
        c("PutativeSubthreshold"),
        c("Control"),
        c("PutativeControl", "Control"),
        c("PutativeControl"),
        c("PutativeControl", "Control"),
        c("PutativeControl", "Control")
    )
    aim2A_results <- lapply(1:length(vars1), function(i_what) {
        var1 <- vars1[[i_what]]
        var2 <- vars2[[i_what]]
        phenoS <- phenoS2
        phenoS$case <- NA
        for(i_v in 1:length(var1)) {
            phenoS$case[phenoS[, "group2018"] == var1[i_v]] <- 1
        }
        for(i_v in 1:length(var2)) {
            phenoS$case[phenoS[, "group2018"] == var2[i_v]] <- 0
        }
        formula1 <- as.formula(paste0("case ~ ", prs, " + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5"))
        formula1_no_prs <- as.formula("case ~ maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")
        ## formula1 <- as.formula("case ~ PRS_PGC_2014_SCZ")
        formula2 <- as.formula(paste0("case ~ ", prs, " + sex + PC1 + PC2 + PC3 + PC4 + PC5")        )
        formula2_no_prs <- as.formula("case ~ sex + PC1 + PC2 + PC3 + PC4 + PC5")
        s1 <- tryCatch({
            (glm(data = phenoS, formula1, family = binomial))
        }, warning = function(w) {
            NA
        })
        ##
        ##s1 <- tryCatch(summary(s1))
        s1B <- tryCatch({
            (glm(data = phenoS, formula1_no_prs, family = binomial))
        }, warning = function(w) {
            NA
        })
        ##s1B <- tryCatch(summary(s1B))
        ##
        s2 <- tryCatch({
            (glm(data = phenoS, formula2, family = binomial))
        }, warning = function(w) {
            NA
        })
        ##s2 <- tryCatch(summary(s2))
        ##
        s2B <- tryCatch({
            (glm(data = phenoS, formula2_no_prs, family = binomial))
        }, warning = function(w) {
            NA
        })
        ##s2B <- tryCatch(summary(s2B))
        ##
        s3 <- tryCatch({
            rms::lrm(formula1, data = phenoS)
        }, warning = function(w) {
            NA
        }, error = function(e) {
            NA
        })
        s3_no_prs <- tryCatch({
            rms::lrm(formula1_no_prs, data = phenoS)
        }, warning = function(w) {
            NA
        }, error = function(e) {
            NA
        })
        s4 <- tryCatch({
            rms::lrm(formula2, data = phenoS)
        }, warning = function(w) {
            NA
        }, error = function(e) {
            NA
        })
        s4_no_prs <- tryCatch({
            rms::lrm(formula2_no_prs, data = phenoS)
        }, warning = function(w) {
            NA
        }, error = function(e) {
            NA
        })
        return(
            list(
                with_age = s1,
                without_age = s2,
                with_age_without_prs = s1B,
                without_age_without_prs = s2B,
                lrm_with_age = s3,
                lrm_with_age_no_prs = s3_no_prs,
                lrm_without_age = s4,
                lrm_without_age_no_prs = s4_no_prs,
                phenoS = phenoS
            )
        )
    })
    names(aim2A_results) <- sapply(1:length(vars1), function(i_what) {
        var1 <- paste0(vars1[[i_what]], collapse = "&")
        var2 <- paste0(vars2[[i_what]], collapse = "&")
        return(paste0(var1, "_vs_", var2))
    })
    return(aim2A_results)
}


aim2A_plot_hist <- function(filename) {
    pdf(filename, height = 4, width = 4)
    hist(
        pheno[, "PRS_PGC_2014_SCZ"],
        xlab = "PRS SCZ",
        ylab = "Count",
        breaks = 20,
        main = ""
    )
    dev.off()
    return(NULL)
}

aim2B_plot_hist <- function(filename) {
    pdf(filename, height = 4, width = 4)
    hist(
        pheno[, "PRS_Neale_XXXX_fluid"],
        xlab = "PRS IQ",
        ylab = "Count",
        breaks = 20,
        main = ""
    )
    dev.off()
    return(NULL)
}

aim2A_plot_groups <- function(filename) {
    pdf(filename, height = 5, width = 7)
    par(mar=c(5.1, 4.1, 1, 2.1), mgp=c(3, 2,0))
    x <- phenoS2[, "PRS_PGC_2014_SCZ"]
    ylim <- range(x)
    ylim[2] <- 0
    ## make ylimit normal range + 20%
    a <- quantile(x, probs = c(0.01, 0.99))
    D <- diff(a)
    minorLevel <- 0.2
    ylim <- c(a[1] - 1 * minorLevel * D, a[2] + 5 * minorLevel * D)
    ybottomlevels <- c(
        a[1] - 1 * minorLevel * D
    )
    ytoplevels <- c(
        a[2] + 1 * minorLevel * D,
        a[2] + 2 * minorLevel * D,
        a[2] + 3 * minorLevel * D
    )
    names <- c("SSD\n", "Subthreshold\npsychosis", "Putative\ncontrol", "Definite\ncontrol")
    boxplot(
        x ~ phenoS2[, "group2018"],
        names = names,
        ylab = "PS_SZ",
        outline = FALSE,
        ylim = ylim,
        col = 'white',
        border = cbPalette,
        axes = FALSE,
        xlab = ""
    )
    ## cols
    axis(1, labels = names, at = 1:4)
    axis(2)
    for(i_group in 1:4) {
        group <- groups_to_plot[i_group]
        w <- phenoS2[, "group2018"] == group
        points(
            x = rnorm(sum(w), mean = i_group, sd = 0.10),
            y = x[w],
            col = scales::alpha(cbPalette[i_group], 0.5),
            pch = 16
        )
    }
    ## Now add p-values and comparisons
    make_label <- function(group, age) {
        p <- coefficients(summary(aim2A_results[[group]][[age]]))["PRS_PGC_2014_SCZ", "Pr(>|z|)"]
        p <- formatC(p, format = "e", digits = 2)
        N <- length(summary(aim2A_results[[group]][[age]])[["deviance.resid"]])
        print((paste0("N = ", N, ", p = ", p)))
        return(paste0("p = ", p))
    }
    ## r2 <- aim2A_results[[1]][["lrm_with_age"]]$stats["R2"]
    ## levels <- list(
    ##     list(c(1.5, 3.5), ybottomlevels[1], make_label("Case_SSD&PutativeSubthreshold_vs_PutativeControl&Control", "with_age")),
    ##     list(c(1, 4),  ytoplevels[2], make_label("Case_SSD_vs_Control", "with_age")),
    ##     list(c(2, 3),  ytoplevels[3], make_label("PutativeSubthreshold_vs_PutativeControl", "with_age")),
    ##     list(c(1, 2), ytoplevels[1], make_label("Case_SSD_vs_PutativeSubthreshold", "with_age")),
    ##     list(c(3, 4), ytoplevels[1], make_label("PutativeControl_vs_Control", "without_age"))
    ## )
    levels <- list(
        list(c(1, 3.5),  ytoplevels[1], make_label("Case_SSD_vs_PutativeControl&Control", "with_age")),
        list(c(1, 4),  ytoplevels[2], make_label("Case_SSD_vs_Control", "with_age"))
    )
    ## list(c(2, 4),  ytoplevels[2], make_label("PutativeSubthreshold_vs_Control", "without_age")),
    for(i in 1:length(levels)) {
        between <- levels[[i]][[1]]
        level <- levels[[i]][[2]]
        label <- levels[[i]][[3]]
        text(x = mean(between), y = level, labels = label, pos = 3)
        arrows(x0 = between[1], x1 = between[2], y0 = level, y1 = level, code = 3, angle = 90, length = 0.05)
        ## if rounded, add two more between span, facing up!
        print(between)
        for(i_b in 1:length(between)) {
            b <- between[i_b]
            if (b != round(b)) {
                segments(x0 = floor(b), x1 = b, y0 = level - minorLevel, y1 = level)
                segments(x0 = b, x1 = ceiling(b), y0 = level, y1 = level - minorLevel)
            }
        }
    }
    dev.off()
    return(NULL)
}


aim2A_plot_r2 <- function(filename) {
    ##
    ## 1 = case_ssd vs control
    ## 5 = putative subthreshold vs control
    ## 4 = both cases vs both controls
    comps <- c("Case_SSD_vs_Control", "PutativeSubthreshold_vs_Control", "Case_SSD&PutativeSubthreshold_vs_PutativeControl&Control")
    fancy_comps <- c("Case_SSD vs\nControl", "PutativeSubthreshold vs\nControl", "Case_SSD &\n PutativeSubthreshold vs \nPutativeControl &\nControl")
    r2A <- array(NA, length(comps))
    r2B <- array(NA, length(comps))
    for(i in 1:length(comps)) {
        p <- comps[i]
        if (! is.na(aim2A_results[[p]]$with_age[1])) {
            r2A[i] <- aim2A_results[[p]][["lrm_with_age"]][["stats"]][["R2"]]
            r2B[i] <- aim2A_results[[p]][["lrm_with_age_no_prs"]][["stats"]][["R2"]]
        } else {
            r2A[i] <- aim2A_results[[p]][["lrm_without_age"]][["stats"]][["R2"]]
            r2B[i] <- aim2A_results[[p]][["lrm_without_age_no_prs"]][["stats"]][["R2"]]
        }
    }
    ##
    ## make barplot here
    ##
    for(suffix in c("pdf", "png")) {
        image_open(filename = filename, height = 4, width = 6, suffix = suffix, remove_suffix = TRUE)
        par(mar=c(5.1, 4.1, 1, 2.1), mgp=c(3, 1,0))
        height <- r2A - r2B
        ##    names(height) <- fancy_comps
        barplot(
            height = height,
            col = cbPalette[4 + 1:length(height)],
            border = NA,
            ylab = "Marginal Nagelkerke R square",
            space = 0.5,
            xlim = c(0.5, 4.5),
            ylim = c(0, 0.09)
        )
        ## add text
        text(
            x = c(1, 2.5, 4),
            y = height,
            labels = round(height, 3),
            pos = 3
        )
        mtext(text = fancy_comps, side = 1, at = c(1, 2.5, 4), padj = 1)
        dev.off()
    }
}


get_aim2B_results <- function(phenotypes) {

    aim2B_results <- lapply(1:length(phenotypes), function(i_what) {

        phenotype <- phenotypes[i_what]
        print(phenotype)
        f <- function(phenotype, prs) {
            return(as.formula(paste0(phenotype, " ~ ", prs, " + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")))
        }
        f2 <- function(phenotype) {
            return(as.formula(paste0(phenotype, " ~ maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5")))
        }
        formula1 <- f(phenotype, "PRS_PGC_2014_SCZ")
        formula2 <- f(phenotype, "PRS_Neale_XXXX_fluid")
        formula1_noPRS <- f2(phenotype)
        formula2_noPRS <- f2(phenotype)
        if (length(grep("binary", phenotype)) > 0) {
            s1 <- glm(data = pheno, formula1, family = binomial)
            s2 <- glm(data = pheno, formula2, family = binomial)
            prs_scz_prediction <- stats::predict(object = s1, newdata = pheno[, rownames(coefficients(summary(s1)))[-1]])
            prs_iq_prediction <- stats::predict.glm(object = s2, newdata = pheno[, rownames(coefficients(summary(s2)))[-1]])
            s1_noPRS <- glm(data = pheno, formula1_noPRS, family = binomial)
            s2_noPRS <- glm(data = pheno, formula2_noPRS, family = binomial)
        } else {
            s1 <- lm(data = pheno, formula1)
            s2 <- lm(data = pheno, formula2)
            ## predict?
            prs_scz_prediction <- stats::predict.lm(object = s1, newdata = pheno[, rownames(coefficients(summary(s1)))[-1]])
            prs_iq_prediction <- stats::predict.lm(object = s2, newdata = pheno[, rownames(coefficients(summary(s2)))[-1]])
            ## also, fit model without PRS
            s1_noPRS <- lm(data = pheno, formula1_noPRS)
            s2_noPRS <- lm(data = pheno, formula2_noPRS)
        }
        return(
            list(
                prs_scz = summary(s1),
                prs_iq = summary(s2),
                prs_scz_prediction = prs_scz_prediction,
                prs_iq_prediction = prs_iq_prediction,
                scz_no_PRS = summary(s1_noPRS),
                iq_no_PRS = summary(s2_noPRS)
            )
        )
    })

    names(aim2B_results) <- phenotypes
    return(aim2B_results)
}

aim2B_plot_groups <- function(filename, what = "pdf", ylimAdjust = 0.2, x.intersp = 0, border = NA, plot_type = "scz") {
    ##if (length(phenotypes) != 4) {
    ##    stop("this plotting assumes a length of 4!")
    ## }
    if (what == "pdf") {
        pdf(filename, height = 4 * 2, width = 4 * 2)
    } else if (what == "png") {
        png(filename, height = 4 * 2, width = 4 * 2, units = "in", res = 300)
    }
    cbo <- 1 ## cbPalette offset, if desired
    par(mfcol = c(2, 2)) ## assumes length(phenotypes)))
    for(i_what in 1:length(phenotypes)) {
        phenotype <- phenotypes[i_what]
        prettyPhenotype <- names(phenotypes)[i_what]
        ## col_start <- cbo + i_what
        col_start <- cbo
        main <- function(prsA, prs2) {
            prs1 <- paste0("prs_", prsA)
            no_prs <- paste0(prsA, "_no_PRS")
            x <- aim2B_results[[phenotypes[i_what]]][[prs1]]
            x_noPRS <- aim2B_results[[phenotypes[i_what]]][[no_prs]]
            beta <- formatC(coefficients(x)[prs2, "Estimate"], digits = 2)
            p <- coefficients(x)[prs2, 4] ## argh different names for p-value
            p <- formatC(p, format = "e", digits = 2)
            N1 <- length(x[["deviance.resid"]])
            N2 <- length(x[["residuals"]])
            N <- max(N1, N2)
            ## also, get r2 difference, specifically linear phenotypes
            r2A <- x$r.squared - x_noPRS$r.squared
            r2B <-
                (1 - x$deviance / x$null.deviance) -
                (1 - x_noPRS$deviance / x_noPRS$null.deviance)
            r2 <- c(r2A, r2B)
            if (length(r2) != 1) {
                stop("bad asumptions!")
            }
            ##return(paste0("N = ", N, "\nr2 = ", round(r2, 4), "\np = ", p))
            return(
                c(
                    paste0("N = ", N),
                    paste0("r2 = ", round(r2, 4)),
                    paste0("p = ", p)
                )
            )
        }
        if ((plot_type == "both") | (plot_type == "scz")) {
            main_scz <- main("scz", "PRS_PGC_2014_SCZ")
        }
        if ((plot_type == "both") | (plot_type == "iq")) {
            main_iq <- main("iq", "PRS_Neale_XXXX_fluid")
        }
        if (length(grep("binary", phenotype)) > 0) {
            f <- function(prs1, xlab, main) {
                xlim <- range(pheno[, prs1], na.rm = TRUE)
                ylim <- c(0.42, 2.58)
                ylim[2] <- ylim[2] + diff(ylim) * ylimAdjust
                boxplot(
                    pheno[, prs1] ~ pheno[, phenotype],
                    xlab = xlab,
                    ylab = prettyPhenotype,
                    main = main,
                    horizontal = TRUE,
                    col = 'white',
                    border = cbPalette[col_start + 1 + 0:1],
                    axes = FALSE,
                    pars = list(
                        xlim = ylim,
                        ylim = xlim
                    )
                )
                usr <- par( "usr" )
                ## add points
                for(i_group in 0:1) {
                    w <- pheno[, phenotype] == i_group
                    w[is.na(w)] <- FALSE
                    y <- rnorm(sum(w), mean = i_group, sd = 0.10)
                    points(
                        x = pheno[w, prs1],
                        y = 1 + y,
                        col = scales::alpha(cbPalette[col_start + 1 + i_group], 0.5),
                        pch = 16
                    )
                }
                axis(1)
                axis(2, at = c(1, 2), labels = c("No", "Yes"))
            }
            if ((plot_type == "both") | (plot_type == "scz")) {
                f("PRS_PGC_2014_SCZ", "PS SZ", "")
                legend("topleft", main_scz, x.intersp = x.intersp, bty = "n")
            }
            if ((plot_type == "both") | (plot_type == "iq")) {
                f("PRS_Neale_XXXX_fluid", "PS IQ", "")
                legend("topleft", main_iq, x.intersp = x.intersp, bty = "n")
            }
            #addText(main_iq, usr, 0.10)
        } else {
            ## quantitative phenotype
            if ((plot_type == "both") | (plot_type == "scz")) {
                x <- pheno[, "PRS_PGC_2014_SCZ"]
                y <- pheno[, phenotype]
                ylim <- range(y, na.rm = TRUE)
                ylim[2] <- ylim[2] + diff(ylim) * ylimAdjust
                plot(x = x, y = y, xlab = "PS SZ", ylab = prettyPhenotype, main = "", col = cbPalette[col_start], axes = FALSE, ylim = ylim)
                axis(1)
                axis(2)
                legend("topleft", main_scz, x.intersp = x.intersp, bty = "n")
                ##            main_scz)
                addTrend(x, y)
                usr <- par( "usr" )
            }
            if ((plot_type == "both") | (plot_type == "iq")) {
                ##
                x <- pheno[, "PRS_Neale_XXXX_fluid"]
                y <- pheno[, phenotype]
                ylim <- range(y, na.rm = TRUE)
                ylim[2] <- ylim[2] + diff(ylim) * ylimAdjust
                plot(x = x, y = y, xlab = "PS IQ", ylab = prettyPhenotype, main = "", col = cbPalette[col_start], axes =FALSE, ylim = ylim)
                axis(1)
                axis(2)
                addTrend(x, y)
                legend("topleft", main_iq, x.intersp = x.intersp, bty = "n")
            }
            ##usr <- par( "usr" )
            ##addText(main_iq, usr)
            ## add simple trendline
            ## add trendline ( WHY did I do all that. just need coefficient, right?
        }
    }
    dev.off()
}

addText <- function(main, usr, v = 0.05) {
    d1 <- v * (usr[2] - usr[1])
    d2 <- v * (usr[4] - usr[3])
    text( usr[2] - d1, usr[4] - d2, main,    adj = c( 1, 1 ))
}
addTrend <- function(x, y) {
    c <- coefficients(summary(lm(y ~ x)))
    abline(coef = c(c[1, 1], c[2, 1]))
}


sum_who <- function(group, group_names) {
    x <- NULL
    for(g in group_names) {
        x <- c(x, sum(group == g))
        names(x)[length(x)] <- g
    }
    return(x)
}


effect_of_cutoff <- function(
    prs,
    groups,
    cutoffs,
    group_names,
    control_groupname,
    altPrev,
    flip_sign = FALSE
) {
    ## make fancy table
    to_out <- array(NA, c(length(cutoffs), length(group_names) + 1))
    if (flip_sign) {
        signs <- c("<", ">=")
    } else {
        signs <- c(">", "<=")
    }
    rownames(to_out) <- c(
        paste0("PRS ", signs[1], " ", cutoffs * 100)
    )
    ##        paste0("PRS ", signs[2], " ", 100 * prs_lower)
    colnames(to_out) <- c(group_names, "all_controls")
    ##
    results_to_return <- as.list(1:length(cutoffs))
    for(i_cutoff in 1:length(cutoffs)) {
        ##
        test_positive <- (!flip_sign) == (quantile(prs, probs = cutoffs[i_cutoff]) <= prs)
        test_negative <- (!flip_sign) == (prs < quantile(prs, probs = cutoffs[i_cutoff]))
        results <- sapply(group_names, function(w) {
            ##
            truth_positive <- array(FALSE, length(prs))
            truth_negative <- array(FALSE, length(prs))
            truth_positive[groups %in% w] <- TRUE
            truth_negative[groups %in% control_groupname] <- TRUE
            broken <- sum(
            (truth_positive & truth_negative) |
            (test_positive & test_negative)
            )
            ##
            tp <- sum(test_positive & truth_positive)
            fp <- sum(test_positive & truth_negative)
            fn <- sum(test_negative & truth_positive)
            tn <- sum(test_negative & truth_negative)
            ##
            or <- (tp / fp) / (fn / tn)
            de <- tp
            he <- fp
            dn <- fn
            hn <- tn
            or2 <- (de / he) / (dn / hn)
            ppv <- tp / (tp + fp)
            ##
            sen <- (tp) / (tp + fn)
            spe <- (tn) / (tn + fp)
            prev <- (tp + fn) / (tp + fn + fp + tn)
            ppv2 <- (sen * prev) / ((sen * prev) + (1 - spe) * (1 - prev))
            ## altPrev <- 0.01
            ppv_altPrev <- (sen * altPrev) / ((sen * altPrev) + (1 - spe) * (1 - altPrev))
            ##
            ## OR CI
            log_or <- log(or)
            log_or_se <- sqrt(1 / tp + 1 / fp + 1 / fn + 1 / tn)
            q <- log_or / log_or_se
            p_or <- 2 * pnorm(q = abs(q), lower.tail = FALSE)
            if (p_or > 1) {
                stop("bad number!")
            }
            ## p-value?
            or_lower <- exp(log(or) - 1.96 * log_or_se)
            or_upper <- exp(log(or) + 1.96 * log_or_se)
            ## PPV CI
            ppv_se <- sqrt((ppv * (1 -ppv) / (tp + fp)))
            ppv_lower <- ppv - 1.96 * ppv_se
            ppv_upper <- ppv + 1.96 * ppv_se
            ##
            to_return <- c(
                tp = tp, fp = fp, fn = fn, tn = tn, or = or,
                or_lower = or_lower,
                or_upper = or_upper,
                p_or = p_or,
                ppv_lower = ppv_lower,
                ppv_upper = ppv_upper,
                ppv = ppv,
                broken = broken,
                ppv_altPrev = ppv_altPrev
            )
            return(to_return)
        })
        ##
        results <- cbind(results, all_controls = rowSums(results[, control_groupname, drop = FALSE]))
        or <- results["or", ]
        or_lower <- results["or_lower", ]
        or_upper <- results["or_upper", ]
        p_or <- results["p_or", ]
        ppv_lower <- results["ppv_lower", ]
        ppv_upper <- results["ppv_upper", ]
        ppv_altPrev <- results["ppv_altPrev", ]
        ppv <- results["ppv", ]
        p <- results["tp", ]
        n <- results["fn", ]
        ## now merge
        ##to_out[i_cutoff, ] <- paste0(round(or, 2), ", (N = ", high, ", ", round(100 * high / sum(high), 1), "%)")
        print(p_or)
        to_out[i_cutoff, ] <- paste0(
            "OR = ", round(or, 2), " ",
            "[", round(or_lower, 2), ", ", round(or_upper, 2), "] p = ", formatC(p_or, format = "e", digits = 2), ", ",
            "PPV = ", round(ppv, 3), ", ",
            "[", round(ppv_lower, 3), ", ", round(ppv_upper, 3), "] , ",
            "N+ = ", p, ", ",
            "N- = ", n, ", ",
            "PPV* = ", round(ppv_altPrev, 3)
        )
        results_to_return[[i_cutoff]] <- results
        ## "(N = ", high, ", ", round(100 * high / sum(high), 1), "%)")
    }
    ##        to_out[nrow(to_out), ] <- paste0("(N = ", low, ", ", round(100 * low / sum(low), 1), "%)")
    ## re-name
    return(
        list(
            to_out = to_out,
            results = results_to_return
        )
    )
}


look_at_scz_quantile <- function(aim2X_table_filename, output_plot = TRUE) {

    ##
    ## here, add vs estimated prevalence, new ppv (ppv*)
    altPrev <- 0.01
    groups <- phenoS2[, "group2018"]
    cutoffs <- c(0.90, 0.75, 0.50)
    control_groupname <- c("PutativeControl", "Control")
    out <- effect_of_cutoff(
        prs = phenoS2[, "PRS_PGC_2014_SCZ"],
        groups = groups,
        cutoffs = cutoffs,
        group_names = c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control"),
        control_groupname = control_groupname,
        flip_sign = FALSE,
        altPrev = altPrev
    )
    to_out <- out$to_out
    ## add single mer

    colnames(to_out) <- c(
        "Schizophrenia spectrum diagnosis",
        "Subthreshold psychosis",
        "Putative control",
        "Control",
        "all_controls"
    )
    write.table(
        to_out,
        file = aim2X_tableSCZ_filename,
        row.names = TRUE,
        col.names = TRUE,
        sep = "\t",
        quote = TRUE
    )
    ## can I add OR CI?
    ## can I add PPV CI?
    ## se(logOR) = sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    ## OR +/- 1.96 * exp(sqrt(1 / a + 1 / b + 1 / c + 1 / d))

    209 / (209 + 216 + 382)
    158 / (158 + 216 + 382)

    ## break down groups by PRS_SCZ?
    cutoffs2 <- (1 - cutoffs)[length(cutoffs):1]
    out2 <- effect_of_cutoff(
        prs = phenoS2[, "PRS_PGC_2014_SCZ"],
        groups = groups,
        cutoffs = cutoffs2,
        group_names = c("Case_SSD", "PutativeSubthreshold", "PutativeControl", "Control"),
        control_groupname = control_groupname,
        altPrev = altPrev,
        flip_sign = TRUE
    )

    to_out <- rbind(out$to_out, out2$to_out)
    colnames(to_out) <- c(
        "Schizophrenia spectrum diagnosis",
        "Subthreshold psychosis",
        "Putative control",
        "Control",
        "all_controls"
    )
    write.table(
        to_out,
        file = aim2X_tableSCZ_filename,
        row.names = TRUE,
        col.names = TRUE,
        sep = ",",
        quote = TRUE
    )


    case_col <- "Case_SSD"
    filename <- file.path(results_dir, "cutoffPSile_scz.png")
    flip_sign <- FALSE
    main <- "Schizophrenia"
    make_PSile_plot(
        out = out,
        out2 = out2,
        filename = filename,
        groups = groups,
        cutoffs = cutoffs,
        cutoffs2 = cutoffs2,
        altPrev = altPrev,
        case_col = case_col,
        control_groupname = control_groupname,
        flip_sign = flip_sign,
        main = main,
        output_plot = output_plot
    )



}

look_at_iq_quantile <- function(aim2X_tableIQ_filename, output_plot = TRUE) {

    groups <- c("ID", "Not ID")[as.integer(phenoS2[, "FSIQ_First"] > 70) + 1]
    prs <- phenoS2[, "PRS_Neale_XXXX_fluid"]
    to_remove <- is.na(groups)
    groups <- groups[!to_remove]
    prs <- prs[!to_remove]
    altPrev <- 0.025

    cutoffs <- c(0.10, 0.25, 0.50) ## these are the upper
    out <- effect_of_cutoff(
        prs = prs,
        groups = groups,
        cutoffs = cutoffs,
        group_names = c("ID", "Not ID"),
        control_groupname = "Not ID",
        flip_sign = TRUE,
        altPrev = altPrev
    )
    to_out <- out$to_out

    ##
    290 / (290 + 411) ## 41%
    ## rest
    44 / (44 + 27) ## 62%
    (290 - 44) / (411 - 27 + 290 - 44) ## 39%

    cutoffs2 <- sort(1 - cutoffs)
    out2 <- effect_of_cutoff(
        prs = prs,
        groups = groups,
        cutoffs = cutoffs2,
        group_names = c("ID", "Not ID"),
        control_groupname = "Not ID",
        flip_sign = FALSE,
        altPrev = altPrev
    )

    to_out <- rbind(out$to_out, out2$to_out)
    write.table(
        to_out,
        file = aim2X_tableIQ_filename,
        row.names = TRUE,
        col.names = TRUE,
        sep = ",",
        quote = TRUE
    )



    main <- "ID (IQ < 70)"
    make_PSile_plot(
        out = out,
        out2 = out2,
        filename = file.path(results_dir, "cutoffPSile.png"),
        groups = groups,
        cutoffs = cutoffs,
        cutoffs2 = cutoffs2,
        altPrev = altPrev,
        case_col = "ID",
        control_groupname = "Not ID",
        flip_sign = TRUE,
        main = main,
        output_plot = output_plot
    )


}





if (1 == 0) {

    ## 22q
    var <- sample(c("22q", "3q"), 1000, replace = TRUE)
    ps <- rnorm(n = 1000)
    iq <-
        70 * (var == "22q") +
        80 * (var == "3q") +
        10 * ps +
        30 * rnorm(n = 1000)
    id <- as.integer((iq < 70))
    ## only using 22q
    data <- data.frame(iq = iq, id = id, ps = ps, var = var)
    data2 <- data[data[, "var"] == "22q", ]
    r <- glm(id ~ ps, data = data, family = "binomial")
    r <- glm(id ~ ps, data = data2, family = "binomial")

    a <- head(ps)
    cbind(
        -0.12437 - 0.65781 * a,
        head(predict(r, data.frame(ps)))
    ) ## OK
    b <- head(predict(r, data.frame(ps)))
    exp(b) / (1 + exp(b)) ## actual predictions

    ##


    x1 <- runif(100)
    x2 <- x1 + runif(100)
    e <- runif(100)
    y <- x2 + e
    z <- y > 1.50

    ## ROC curve y ~ x2
    ## ROC curve y ~ x2 + x1
    summary(lm(y ~ x2))
    summary(lm(y ~ x2 + x1))
    summary(lm(y ~ x1))



}

plotPointsForPrevalencePercentiles <- function(
    x,
    m,
    m2,
    col,
    control_cols,
    column = "Case_SSD",
    w = 0.1
) {
    y <- 100 * m2[, column]
    points(x = x, y = y, type = "o", col = col, pch = 16) ## case
    ##ci <- 100 * (1.96 * sqrt(m2[, column] * (1 - m2[, column]) / (m[, column] + m[, 3] + m[, 4])))
    ci <- 100 * (1.96 * sqrt(m2[, column] * (1 - m2[, column]) / (rowSums(m[, c(column, control_cols)]))))
    arrows(x0 = x, x1 = x, y0 = y - ci, y1 = y + ci, length = 0, col = col) ## bars
    arrows(x0 = x - w, x1 = x + w, y0 = y + ci, y1 = y + ci, length = 0, col = col) ## bars
    arrows(x0 = x - w, x1 = x + w, y0 = y - ci, y1 = y - ci, length = 0, col = col) ## bars
}


plotPrevalenceByPercentiles <- function(
    file,
    prs,
    phenoS2,
    control_cols,
    plot_groups,
    fancy_plot_groups,
    category_column = "group2018",
    ylim1 = c(0, 35),
    ylim2 = c(0, 2),
    quantile_probs = c(0, 0.25, 0.5, 0.75, 1.0),
    altPrev = 0.01
) {
    ## labels <- c("[0-25)", "[25-50)", "[50-75)", "[75-100]")
    a <- 100 * quantile_probs
    labels <- paste0("[", a[-length(a)], "-", a[-1], ")")
    labels[length(labels)] <- gsub(")", "]", labels[length(labels)])
    ##
    q <- quantile(prs, probs = quantile_probs)
    m <- table(cut(prs, q), phenoS2[, category_column])
    m2 <- round(m / rowSums(m), 2)
    m2ORI <- m2
    ##
    for(iAltPrev in 1:2) {
        if (iAltPrev == 2) {
            file <- gsub(".png", ".altPrev.png", file)
        }
        ylim <- list(ylim1, ylim2)[[iAltPrev]]
        png(file, height = 6, width = 6, res = 300, units = "in")
        cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
        xlim <- c(-0.1, 1.1)
        ##    ylim <- c(0, 0.7)
        ##x <- seq(0, 1, length.out = 4)
        mids <- (quantile_probs[-1] + quantile_probs[-length(quantile_probs)]) / 2
        x <- mids
        plot(x = 0, y = 0, xlim = xlim, ylim = ylim, xlab = "Quartile", ylab = "Frequency (%)", col = "white", axes = FALSE)
        axis(1, at = x, labels = labels)
        axis(2)
        ## case
        for(i_plot_group in 1:length(plot_groups)) {
            ##
            m2 <- m2ORI
            plot_group <- plot_groups[i_plot_group]
            ##
            if (iAltPrev == 2) {
                oriPrev <- sum(m[, plot_group]) / sum(m[, c(plot_group, control_cols)])
                ## now, make new m
                m2[, plot_group] <- m2[, plot_group] / oriPrev * altPrev
            }
            ## here, change!
            plotPointsForPrevalencePercentiles(x, m, m2, col = cbPalette[1 + i_plot_group], column = plot_group, w = 0.1, control_cols = control_cols)
        }
        ##    points(x = x, y = m2[, 3] + m2[, 4], type = "l", col = cbPalette[4]) ## case
        legend("bottomright", fancy_plot_groups, col = cbPalette[1 + 1:length(plot_groups)], lwd = 2)
        ## add CI
        dev.off()
    }
}


plot_percentiles_vs_prevalence <- function() {

    quantile_probs <- c(0, 0.25, 0.5, 0.75, 1)
    quantile_probs <- seq(0, 1, length.out = 6)
    ## has PPV as well
    file <- file.path(results_dir, "scz.prev_by_prs_tile.png")
    prs <- phenoS2[, "PRS_PGC_2014_SCZ"]
    plot_groups <- c("Case_SSD") #, "PutativeSubthreshold")
    control_cols <- c("PutativeControl", "Control")
    category_column <- "group2018"
    fancy_plot_groups <- c("Schizophrenia") #, "Subthreshold psychosis")
    ylim1 <- c(0, 35)
    ylim2 <- c(0, 2)
    altPrev <- 0.01
    plotPrevalenceByPercentiles(
        file = file,
        prs = prs,
        phenoS2 = phenoS2,
        category_column = category_column,
        control_cols = control_cols,
        plot_groups = plot_groups,
        fancy_plot_groups = fancy_plot_groups,
        ylim1 = ylim1,
        ylim2 = ylim2,
        quantile_probs = quantile_probs,
        altPrev = altPrev
    )
    ## make another plot afterwards, use rr vs average, with percentile, make new values?
    ## ID
    file <- file.path(results_dir, "id.prev_by_prs_tile.png")
    prs <- phenoS2[, "PRS_Neale_XXXX_fluid"]
    phenoS3 <- phenoS2
    phenoS3$ID <- c("Not ID", "ID")[as.integer(phenoS2[, "FSIQ_First"] < 70) + 1]
    plot_groups <- c("ID")
    control_cols <- c("Not ID")
    category_column <- "ID"
    fancy_plot_groups <- "ID"
    ylim1 <- c(0, 60)
    ylim2 <- c(0, 6)
    altPrev <- 0.025
    plotPrevalenceByPercentiles(
        file = file,
        prs = prs,
        phenoS2 = phenoS3,
        category_column = category_column,
        control_cols = control_cols,
        plot_groups = plot_groups,
        fancy_plot_groups = fancy_plot_groups,
        ylim1 = ylim1,
        ylim2 = ylim2,
        quantile_probs = quantile_probs,
        altPrev = altPrev
    )
    ##

    ## PPV (prevalence) for other groups

}



old_code <- function() {


    tp <- sum(groups == "ID" & prs < quantile(prs, 0.10))
    fn <- sum(groups == "ID" & prs > quantile(prs, 0.10))
    fp <- sum(groups != "ID" & prs < quantile(prs, 0.10))
    tn <- sum(groups != "ID" & prs > quantile(prs, 0.10))
    ppv <- tp / (tp + fp)
    npv <- tn / (tn + fn)
    sen <- tp / (tp + fn)
    spe <- tn / (tn + fp)
    c(ppv = ppv, npv = npv, sen = sen, spe = spe)


## load up into RAM?
pdf_file <- file.path(results_dir, "hist.pdf")
pdf(pdf_file, height = 8, 12)
par(mfrow = c(2, 3))
for(phenotype in names(paths_to_phenotypes)) {
    all_score_file <- file.path(results_dir, paste0(phenotype, "_iBBC"), "PRSice.all.score")
    all_score <- read.table(all_score_file, header = TRUE)
    ## just choose this one for now
    hist(all_score[, "X0.100000"], main = pretty_phenotype_names[phenotype], xlab = "Value")
}
dev.off()
## show distributions

## specifically regress SCZ 0.1 PRS against SCZ, and intelligence 0.1 against IQ
pdf_file <- file.path(results_dir, "test.pdf")
pdf(pdf_file, height = 4, width = 8)
par(mfrow = c(1, 2))
## schizophrenia
model <- glm(case ~ PRS, family = binomial(link = 'logit'),data=phenoS)
summary(model)
x <- sort(phenoS$PRS[phenoS$case == 1])
plot(x = seq(0, 1, length.out = length(x)), y = x, col = cbPalette[2], main = "PGC 2014 SCZ", ylab = "PRS", xlab = "Rank")
points(
    seq(0, 1, length.out = sum(phenoS$case == 0, na.rm = TRUE)),
    y = sort(phenoS$PRS[phenoS$case == 0]), col = cbPalette[3]
)
legend("topright", c("Case_SSD", "Control"), col = cbPalette[2:3], lwd = 2)
## intelligence
phenoS <- add_prs_to_pheno(phenotype = "sniekers_2017_intelligence")
summary(lm(FSIQ_First ~ PRS, data = phenoS))
## scatterplot
r2 <- cor(phenoS[, "PRS"], phenoS[, "FSIQ_First"], use = "pairwise.complete") ** 2
plot(x = phenoS[, "PRS"], y = phenoS[, "FSIQ_First"], xlab = "PRS (Sniekers 2017)", ylab = "FSIQ_First", main = paste0("Sniekers 2017 Intelligence\nr2 = ", round(r2, 3)))
##
dev.off()


quit()

## z-score, p-value
## linear regression
set.seed(10)
x <- runif(1000)
y <- 0.2 * x + runif(1000)
summary(lm(y ~ x))

##
davies <- data.table::fread(davies_2018_gcf, data.table = FALSE)
neale <- data.table::fread(paste0("gunzip -c ", file.path(external_dir, "fluid_intelligence.20016.assoc.tsv.gz")), data.table = FALSE)

## compare p-values for overlapping SNPs?
both <- intersect(neale[, "rsid"], davies[, "MarkerName"])
nealeI <- neale[match(both, neale[, "rsid"]), ]
daviesI <- davies[match(both, davies[, "MarkerName"]), ]

## correlation between p-values
cor(nealeI[, "pval"], daviesI[, "P"])

x <- -log10(nealeI[, "pval"])
y <- -log10(daviesI[, "P"])
w <- (2 < x) | (2 < y)

png("~/fluid.png")
plot(x[w], y[w], xlab = "Neale", ylab = "Davies")
dev.off()

}




for_jacob_plot_fraction_of_controls_and_ps_scz <- function() {

    ## average among groups
pheno <- read.csv(pheno_file_with_prs)
    phenoS2 <- make_pheno_with_ordered_group2018(pheno)
    x1 <- mean(phenoS2[phenoS2[, "group2018"] == "Case_SSD", "PRS_PGC_2014_SCZ"])
    x2 <- mean(phenoS2[phenoS2[, "group2018"] == "Control", "PRS_PGC_2014_SCZ"])
    x3 <- mean(phenoS2[phenoS2[, "group2018"] == "PutativeSubthreshold", "PRS_PGC_2014_SCZ"])
    x4 <- mean(phenoS2[phenoS2[, "group2018"] == "PutativeControl", "PRS_PGC_2014_SCZ"])
    ##
    xx <- seq(0, 1, length.out = 100)
    vals <- sapply(xx, function(x) {
        (1 - x) * x1 + (x) * x2
    })
    ##
    ## simulation for Jacob
    ## simulate group
    ci_f <- function(y) {
        m <- mean(y)
        ci <- sd(y) / sqrt(length(y))
        ci <- c(m - 1.96 * ci, m + 1.96 * ci)
        return(ci)
    }
    find_match <- function(x, xx, vals) {
        return(xx[which.min(abs(x - vals))])
    }
    ##
    pdf(file.path(results_dir, "cont.frac.ps_scz.pdf"), height = 6, width = 6)
    xlim <- c(-0.1, 1.1)
    ylim <- c(-0.3, 0.4)
    plot(xx, vals, xlab = "Fraction of controls", type = "l", ylab = "PS_SZ", ylim = ylim, xlim = xlim, axes = FALSE)
    ## add end bits
    points(x = c(head(xx, 1), tail(xx, 1)), y = c(head(vals, 1), tail(vals, 1)), type = "o", cex = 2)
    axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), labels = c("SSD value", 0.2, 0.4, 0.6, 0.8, "Ctrl value"))
    axis(2)
    ##
    ## plot where subthrehsold cases would be
    ##text(x = 0, y = x1, "Schizophrenia\nvalue", col = "red")
    ##text(x = 1, y = x2, "Control\nvalue", col = "red")
    ## subthreshold value
    add_to_plot <- function(phenotype, print_pheno, xx, vals) {
        av_for_pheno <- mean(phenoS2[phenoS2[, "group2018"] == phenotype, "PRS_PGC_2014_SCZ"])
        x <- xx[which.min(abs(av_for_pheno - vals))]
        ci <- ci_f(phenoS2[phenoS2[, "group2018"] == phenotype, "PRS_PGC_2014_SCZ"])
        arrows(x0 = x, x1 = x, y0 = ci[1], y1 = ci[2], col = "purple", code = 3, angle = 90)
        ##
        spot_for_pheno <- find_match(av_for_pheno, xx, vals)
        ## lower ci point
        a <- round(find_match(ci[1], xx, vals), 2)
        text(x = spot_for_pheno, y = ci[1], a, pos = 1)
        ## upper ci point
        a <- round(find_match(ci[2], xx, vals), 2)
        text(x = spot_for_pheno, y = ci[2], a, pos = 3)
        ##
        a <- round(find_match(av_for_pheno, xx, vals), 2)
        text(x = spot_for_pheno, y = av_for_pheno, a, pos = 2)
        ## print_pheno <- paste0(print_pheno, "\n% = ", round(100 * spot_for_pheno, 1))
        print_pheno <- paste0(print_pheno) ## , "\n% = ", round(100 * spot_for_pheno, 1))
        ## move away a bit
        wx <- diff(xlim) * 0.2
        wy <- diff(ylim) * 0.1
        text(x = spot_for_pheno + wx, y = av_for_pheno + wy, print_pheno, col = "red")
        ## grey line for joining text to red dot
        lines(x = c(spot_for_pheno, spot_for_pheno + wx), y = c(av_for_pheno, av_for_pheno + wy), col = "grey")
        ## red dor
        points(x = spot_for_pheno, y = av_for_pheno, col = "red", cex = 1.5, pch = 16)
        ## grey line down to bottom
    }
    ##
    add_to_plot(phenotype = "PutativeSubthreshold", print_pheno = "Subthreshold\npsychosis", xx, vals)
    add_to_plot(phenotype = "PutativeControl", print_pheno = "Putative\nControl", xx, vals)
    ## grey lines
    lines(x = c(0, 0), y = c(-1 + ylim[1], head(vals, 1)), col = "grey")
    lines(x = c(1, 1), y = c(-1 + ylim[1], tail(vals, 1)), col = "grey")
    dev.off()
    ##

}


