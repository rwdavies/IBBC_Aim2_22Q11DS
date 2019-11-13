mediation <- function() {

    ## among subthrehsold cases + putative controls + sex + PC1-5
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
    qsp[, "SIPS1"] <- log(qsp[, "SIPS1"] + 0.5)
    qsp[, "SIPS1"] <- (qsp[, "SIPS1"] - mean(qsp[, "SIPS1"])) / sd(qsp[, "SIPS1"])
    hist(qsp[, "SIPS1"], breaks = 20)
    both <- merge(qsp, phenoS3, by = "genomics_id")

    
    mood <- read.csv("~/IBBC/external/AIMII_mood_2018_01_22.csv")
    both <- mood
    both$case <- NA
    both$case[both[, "group2018"] == "PutativeSubthreshold"] <- 1
    both$case[both[, "group2018"] == "PutativeControl"] <- 0
    both$case[both[, "group2018"] == "Control"] <- 0
    both[both[, "psy_mood"] == "unk", "psy_mood"] <- NA
    ## N = 540
    f <- function(formula1, both, what = "logistic") {
        if (what == "logistic") {
            results <- summary(glm(data = both, as.formula(formula1), family = binomial))
        } else if (what == "linear") {
            results <- summary(lm(data = both, as.formula(formula1)))
        }
        ##print(coefficients(results))
        print(coefficients(results)[2, ])
        ## print(coefficients(results)["PRS_PGC_2014_SCZ", ])        
        print(paste0("N = ", sum(results$df[1:2])))
    }
    ## 1 = mediator and IV
    f("psy_mood ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5", both, "logistic")
    f("psy_mood ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5", both[is.na(both[, "case"]) == FALSE, ], "logistic")    
    ## 2 = DV and IV, absence of mediator
    f("case ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5", both, "logistic")
    ## 3 = DV and mediator, no IV
    f("case ~ psy_mood + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5", both, "logistic")
    ## 4 = DV and IV, presence of mediator
    f("case ~ PRS_PGC_2014_SCZ + psy_mood + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5", both, "logistic")
    ## f("case ~ PRS_PGC_2014_SCZ + maxassessmentage + sex + PC1 + PC2 + PC3 + PC4 + PC5", both, "logistic")


}
