## OK - this works!
## tomorrw - 

R_dir <- "~/proj/22Q11/R/"; results_dir <- "/data/smew1/rdavies/22Qresults/"
source(file.path(R_dir, "functions.R"))
source(file.path(R_dir, "simulate_functions.R"))

## launch jobs here
library("parallel")
N <- 8 * 16
filenamer <- function(i_repeat) {
    return(file.path(results_dir, paste0("learn2.", i_repeat, ".RData")))
}

to_run <- 1:N
to_run <- (4*16 + 1):N
out <- mclapply(to_run, mc.cores = 16, function(i_repeat) {
    file <- tempfile()
    output_file <- filenamer(i_repeat)
    cat(
        'R_dir <- "~/proj/22Q11/R/"; results_dir <- "/data/smew1/rdavies/22Qresults/"',
        'library("VGAM"); library("parallel")',
        'source(file.path(R_dir, "functions.R"))',
        'source(file.path(R_dir, "simulate_functions.R"))',
        paste0('set.seed(', i_repeat, ')'),
        'model_params <- get_validation_model_params()',
        'pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno; ',
        'model_out <- fit_model(pheno, verbose = FALSE, n_ci_iterations = 10, calculate_ci = TRUE, nCores = 1, ci_verbose = TRUE, age_size = model_params$age_size)',
        paste0('save(pheno, model_out, file = "', output_file, '")'),
        sep = "\n",
        file = file
    )
    log_file <- paste0(output_file, ".log")
    system(paste0("R -f  ", file, " &> ", log_file))
})

## am here
## 
## load(file = file.path(results_dir, paste0("learn.", i_repeat, ".RData")))
quit()

model_params <- get_validation_model_params()
vals <- c("K_scz"  ,      "K_sub"       , "age_shape1",   "age_shape2",   "scz_mean_age", "scz_sd_age",   "sub_mean_age", "sub_sd_age")
out <- lapply(1:N, function(i_repeat) {
    ## seems OK so far! if anything, too lose
    file <- filenamer(i_repeat) 
    if (file.exists(file)) {
        load(file = file)
        ## check inside or outside?
        ci_matrix <- model_out$ci_matrix
        if (length(ci_matrix) > 2) {
            contained <- sapply(1:nrow(ci_matrix), function(i_ci) {
                ## val <- model_params[[rownames(ci_matrix)[i_ci]]]
                val <- model_params[[vals[i_ci]]]                
                contained <- (ci_matrix[i_ci, 1] <= val) & (val <= ci_matrix[i_ci, 2])
                return(contained)
            })
        } else {
            contained <- rep(FALSE, 8)
        }
        estimate <- model_out$results$par
        return(list(estimate = estimate, contained = contained, ci_matrix = ci_matrix))
    } else {
        return(NA)
    }
})
out2 <- out[sapply(out, length) > 1]

## get count!
out3 <- sapply(out2, function(x) x$contained)
message("The percent of times the CI is right is:")
out3B <- round(100 * rowSums(out3) / ncol(out3), 1)
out3B

message("The estimate are:")
out4 <- rbind(
    apply(sapply(out2, function(x) x$estimate), 1, mean),
    apply(sapply(out2, function(x) x$estimate), 1, sd)
)
## also, "average ci"?
out4 <- rbind(unlist(model_params[match(colnames(out4), names(model_params))]), out4)
print(ncol(out3))
out4B <- round(rbind(out4, out3B), 3)
print(out4B)
rownames(out4B) <- c("Truth", "Average estimate", "SD estimate", "% times CI contains truth")
write.table(
    out4B,
    file = "~/proj/22Q11/ci.validation.csv",
    row.names = TRUE,
    col.names = TRUE,
    sep = ",",
    quote = FALSE
)

## 
apply(t(sapply(out2, function(x) x$ci_matrix[1, ])), 2, mean) ## average CI


N <- 1e5
age_shape2 <- 4
age <- rbetabinom.ab(n = N, size = age_size, shape1 = age_shape1, shape2 = age_shape2) + 5
## check? quikcly?
for(age_shape2  in c(3, 4, 5)) {
    s <- seq(1.9, 2.1, length.out = 21)
ll <- sapply(s, function(age_shape1) {
    p_age <- dbetabinom.ab(x = round(age - 5), size = age_size, shape1 = age_shape1, shape2 = age_shape2)
    sum(log10(p_age))
})
    print(max(ll))
b <- ll - max(ll)
names(b) <- s
b
}
