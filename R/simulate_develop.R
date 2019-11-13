library("VGAM")
library("rmutil")
library("parallel")
library(gridExtra)
library(lattice)
library(raster)
library(rasterVis)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

R_dir <- get_and_sanitize("R_DIR")
results_dir <- get_and_sanitize("RESULTS_DIR")

if (1 == 1) {

    ## R_dir <- "~/proj/22Q11/R/"; results_dir <- "/data/smew1/rdavies/22Qresults/"; nCores <- 16
    R_dir <- "~/Dropbox/22Q11/R/";    results_dir <- file.path("~/IBBC/", "2018_11_28"); nCores <- 4
    
}


source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "simulate_functions.R"))

set.seed(7958) ## ## GE share volume Nov 27 2018 remove last two, 70,795,829

## test the math here
model_params <- get_default_model_params()
model_params$N <- 100000 ## choose larger N for testing
out <- simulate_full(model_params = model_params, do_checks = TRUE)$pheno

## fit using real data!
file <- file.path(results_dir, "model_params.RData")
if (!file.exists(file)) {
    out <- fit_model_real_data(
        pheno_file = file.path(results_dir, "iBBC_AIMIIdata_14June2018.withPRS.csv")
    )
    model_params <- out$estimated_model_params
    save(model_params, file = file)
} else {
    load(file = file)
}


## simulate an example
pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno

analyze_pheno_for_Aim2A(pheno, group2018_name = "group2018") ## same as for real variable
analyze_pheno_for_Aim2A(pheno, group2018_name = "group2018_rb")

analyze_pheno_for_Aim2AB(pheno, subpheno = "sub", group2018_name = "group2018")
analyze_pheno_for_Aim2AB(pheno, subpheno = "subQ", group2018_name = "group2018")
analyze_pheno_for_Aim2B(pheno)

n_power_reps <- 200
nGrids <- 21
source(file.path(R_dir, "simulate_functions.R"))

## aim2B
power_analysis_aim2B(
    fileprefix = file.path(results_dir, "aim2b.power"),
    n_power_reps = n_power_reps,
    nGrids = nGrids,
    model_params = model_params
)

## for "aim2A", i.e. subthreshold SCZ, do analysis here
power_analysis_aim2A(
    filename = file.path(results_dir, "aim2a.power"),
    n_power_reps = n_power_reps,
    group2018_name = "group2018",
    model_params = model_params,
    prs_colname = "PRS_scz"
)

## flip it as well
power_analysis_aim2A(
    filename = file.path(results_dir, "aim2a.iq.power"),
    n_power_reps = n_power_reps,
    group2018_name = "group2018",
    model_params = model_params,
    prs_colname = "PRS_iq"
)

## do same thing for "random binary" phenotype
power_analysis_aim2A(
    filename = file.path(results_dir, "aim2a.power.rb"),
    n_power_reps = n_power_reps,
    group2018_name = "group2018_rb",
    model_params = model_params
)


## are here!
power_analysis_aim2AB(
    fileprefix = file.path(results_dir, "aim2ab"),
    n_power_reps = 10000,
    group2018_name = "group2018",
    model_params = model_params
)


power_plots_aim2A() ## make single aim2A plot
make_simple_power_matrix() ## make single table with what I want





quit()


system(paste0("rsync -av /data/smew1/rdavies/22Qresults/* rescomp:/well/myers/rwdavies/22Qresults/"))

## 

for(i in 1:6) {
    x <- out_mat[, , i] / c_mat
    levelplot(x, main = dimnames(out_mat)[i])
}

average_PRS_scz <- c(
    Z_scz = mean(G_l1[Z_scz]),
    Z_sub = mean(G_l1[Z_sub]),
    Z_putcont = mean(G_l1[Z_putcont]),
    Z_cont = mean(G_l1[Z_cont])
)
    


## check PRS in subjects
mean(results$data[, "PRS_scz"])
mean(results$data[, "PRS_scz"])
mean(G_l1[Z_sub])
mean(G_l1[Z_putcont])
mean(G_l1[Z_cont])

quit()

## more scratch
pheno <- read.csv("~/IBBC/2018_11_15/iBBC_AIMIIdata_14June2018.withPRS.csv", header = TRUE)

par(mfrow = c(2, 3))
t <- table(pheno[, "group2018"])
for(n in names(t)) {
    w <- pheno[, "group2018"] == n
    plot(sort(pheno[w, "maxassessmentage"]), main = n)
}
