set.seed(4345)

pheno <- read.csv("~/iBBC/external/iBBC_AIMIIdata_14June2018.csv")
table(pheno[, "group2018"])
table(pheno[, "maxassessmentage"], pheno[, "group2018"])[, -1]

N <- nrow(pheno)
get_K <- function(who, var = "group2018") {
    return(sum(pheno[, var] == who, na.rm = TRUE) / sum(is.na(pheno[, var]) == FALSE))
}
#K_case <- get_K("Case_SSD")
#K_control <- get_K("Control")
#K_PutativeSubthreshold <- get_K("PutativeSubthreshold")
#K_PutativeControl <- get_K("PutativeControl")
#K_male <- get_K(1, "sex")
#age <- pheno[, "maxassessmentage"] ## not ideal? but OK for now
K_case <- 0.15 ## lifetime risk
K_subthreshold <- 0.XXX
rm(pheno)

h2_g_SCZ_known <- 0.05 ## explain 5% of variance for liability using SCZ PRS
case_z_thresh <- qnorm(p = 1 - K_case)
PutativeSubthreshold_thresh <- 

G_SCZ_known <- rnorm(n = N, mean = 0, sd = sqrt(h2_g_SCZ_known))
G_SCZ_unknown <- rnorm(n = N, mean = 0, sd = sqrt(1 - h2_g_SCZ_known))
G_SCZ <- G_SCZ_known + G_SCZ_unknown
Z_SCZ <- G_SCZ > case_z_thresh ## by definition
## make subthresold phenotype 50% correlated 
G_subthreshold <- sqrt(0.5) * G_SCZ + rnorm(n = N, mean = 0, sd = sqrt(0.5))
Z_subthreshold <- G_subthreshold > K_PutativeSubthreshold

## assumes sigmoid function fit with
## 10% chance of SCZ at 16
## 90% chance of SCZ at 25
p_SCZ_given_age <- function(age) {
    age1 <- 15
    age2 <- 25
    prob1 <- 0.10
    prob2 <- 0.99
    A1 <- log(1 / prob1 - 1)
    A2 <- log(1 / prob2 - 1)
    b <- (age2 * A1 - age1 * A2) / (A1 - A2)
    A <- -(A2 / (age2 - b))
    return(1 / (1 + exp(-A * (age - b))))
}

## 
group2018 <- array(NA, N)


p_SCZ <- sapply(age, p_SCZ_given_age) ## probability of expressing that phenotype
## cases meet case definition, regardless of age
group2018[Z_SCZ == 1] <- "Case_SSD"
## subthreshold cases have subthreshold phenotype 
[Z_subthreshold == 1 & (age < 25 & SCZ == 0)]
## controls are those with age >25 and no strict phenotype
condition1 <- (Z_SCZ == 0) & (age > 25) ## for now, downsample
group2018[runif(N) < (K_control * N / sum(condition1, na.rm = TRUE)) & (condition1)] <- "Control"
## putative controls are those age<25 without yet phenotype


## write output
z <- array(NA, N)
z[group2018 == "Case_SSD"] <- 1
z[group2018 == "Control"] <- 0

sim_pheno <- cbind(
    group2018 = group2018,
    PRS_PGC_2014_SCZ = G_SCZ_known,
    PC1 = rnorm(N),
    PC2 = rnorm(N),
    PC3 = rnorm(N),
    sex = 1 + as.integer(runif(N) > (1 - K_male))
)

pheno <- write.table(
    sim_pheno,
    file = "~/iBBC/2018_06_18/sim.csv",
    row.names = FALSE
)
