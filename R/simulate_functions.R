ll_f <- function(
    theta,
    pheno,
    age_col = "age",
    age_min = 10,
    age_max = 70,
    mult = 1,
    verbose = FALSE,
    age_size = NA,
    fixed_name = NA,
    fixed_value = NA
) {
    ##
    if (verbose) {
        print(theta)
    }
    if (!is.na(fixed_name)) {
        theta[fixed_name] <- fixed_value
    }
    K_scz <- theta["K_scz"]
    K_sub <- theta["K_sub"]
    ## age_size <- theta["age_size"]
    age_shape1 <- theta["age_shape1"]
    age_shape2 <- theta["age_shape2"]
    scz_mean_age <- theta["scz_mean_age"]
    scz_sd_age <- theta["scz_sd_age"]
    sub_mean_age <- theta["sub_mean_age"]
    sub_sd_age <- theta["sub_sd_age"]
    ## given counts, calculate likelihood
    ## do all at once
    age <- pheno[, age_col]
    norm_age <- (age - age_min) / (age_max - age_min)
    group <- pheno[, "group2018"]
    raw_age <- age - 5
    if (is.na(age_size)) {
        age_size <- max(raw_age)
    }
    ##print(age_size)
    ##stop("WER")        
    if (sum(raw_age < 0) > 0) {
        stop("Bad age")
    }
    p_age <- dbetabinom.ab(x = raw_age, size = age_size, shape1 = age_shape1, shape2 = age_shape2)
    ## probability of PRS given NOTHING
    prs <- pheno[, "PRS_scz"]    
    p_prs <- dnorm(prs, mean = 0, 1)
    ## 
    local_scz_prob <- K_scz * pnorm(q = age, mean = scz_mean_age, sd = scz_sd_age)
    local_sub_prob <- K_sub * pnorm(q = age, mean = sub_mean_age, sd = sub_sd_age)
    ## setup
    local_prob <- array(NA, nrow(pheno))
    ##
    probs <- cbind(
        Case_SSD = local_scz_prob,
        Control = (1 - local_scz_prob),
        PutativeSubthreshold = (1 - local_scz_prob) * local_sub_prob,
        PutativeControl = (1 - local_scz_prob) * (1 - local_sub_prob)
    )
    ## age is irrelevant here?
    ##w1 <- (age >= 25)
    ## and case -> case prob
    w2 <- (group == "Case_SSD")
    local_prob[w2] <- probs[w2, "Case_SSD"]
    ## and control -> control prob
    w2 <- (group == "Control")
    local_prob[w2] <- probs[w2, "Control"]
    ## le 25
    ##w1 <- (age < 25)    
    ## and put sub -> sub
    w2 <- (group == "PutativeSubthreshold")
    local_prob[w2] <- probs[w2, "PutativeSubthreshold"]
    ## and put cont -> cont
    w2 <- (group == "PutativeControl")
    local_prob[w2] <- probs[w2, "PutativeControl"]
    ll <- mult * sum(log(local_prob) + log(p_age))
    if (verbose) {    
        message(ll)
    }
    if (is.na(ll)) {
        ll <- 1e10
    }
    if (ll == Inf) {
        ll <- 1e10
    }
    if (ll == -Inf) {
        ll <- 1e10
    }
    return(ll)
    ## over 25 and not case -> control
    ## w <- (age > 25) & (group == "Case_SSD")
    ## local_prob[w] <- (1 - local_scz_prob)[w]
    ## ##
    ## ## under 25 and case -> case
    ## w <- (age <= 25) & (group == "Case_SSD")
    ## local_prob[w] <- local_scz_prob[w]
    ## ## under 25 and not case, but yes put -> put
    ## w <- (age <= 25) & (group == "PutativeSubthreshold")
    ## local_prob[w] <- ((1 - local_scz_prob) * local_sub_prob)[w]
    ## ## under 25 and not cases and not putative -> control
    ## w <- (age <= 25) & (group == "PutativeSubthreshold")    
    ## local_prob[w] <- ((1 - local_scz_prob) * (1 - local_sub_prob))[w]
    ## ##
    ##
    ## for(i in 1:nrow(pheno)) {
    ##     age <- pheno[i, age_col]
    ##     norm_age <- (age - age_min) / (age_max - age_min)
    ##     group <- pheno[i, "group2018"]
    ##     ## probability of age
    ##     ##        p_age <- dbeta(x = norm_age, shape1 = age_shape1, shape2 = age_shape2)
    ##     p_age <- dbetabinom.ab(x = round(age), size = age_size, shape1 = age_shape1, shape2 = age_shape2)
    ##     ## phenotype
    ##     local_scz_prob <- K_scz * pnorm(q = age, mean = scz_mean_age, sd = scz_sd_age)
    ##     local_sub_prob <- K_sub * pnorm(q = age, mean = sub_mean_age, sd = sub_sd_age)
    ##     ## make all options
    ##     if (age > 25) {
    ##         if (group == "Case_SSD") {
    ##             local_prob <- local_scz_prob
    ##         } else {
    ##             local_prob <- 1 - local_scz_prob
    ##         }
    ##     } else {
    ##         if (group == "Case_SSD") {
    ##             local_prob <- local_scz_prob
    ##         } else if (group == "PutativeSubthreshold") {
    ##             local_prob <- (1 - local_scz_prob) * local_sub_prob
    ##         } else {
    ##             local_prob <- (1 - local_scz_prob) * local_sub_prob
    ##         }
    ##     }
    ##     local_ll <- log(p_age) + log(local_prob)
    ##     if (is.nan(local_ll)) {
    ##         stop(i)
    ##     }
    ##     ll <- ll + local_ll
    ## }
}


fit_model <- function(
    pheno,
    K_scz = 0.3,
    K_sub = 0.3,
    age_shape1 = 2,
    age_shape2 = 5, ## does this not matter?
    scz_mean_age = 20,
    scz_sd_age = 4,
    sub_mean_age = 10,
    sub_sd_age = 4,
    verbose = FALSE,
    nCores = 1,
    n_ci_iterations = 1,
    calculate_ci = FALSE,
    ci_verbose = FALSE,
    age_size = NA
) {
    ## 
    theta <- c(
        K_scz = K_scz,
        K_sub = K_sub,
        age_shape1 = age_shape1,
        age_shape2 = age_shape2,
        scz_mean_age = scz_mean_age,
        scz_sd_age = scz_sd_age,
        sub_mean_age = sub_mean_age,
        sub_sd_age = sub_sd_age
    )
    ## simulate theta within range?
    p <- length(theta)
    k <- 2 * p
    ui <- array(0, c(k, p))
    ci <- c(
        0.05, -0.7, ## K_scz > T, -(T) > K_scz
        0.05, -0.7, ## K_sub
        0.1, -20, ## age shape 1
        0.1, -20, ## age shape 2,
        1, -50, ## scz mean
        0.01, -30, ## scz sd
        1, -50, ## sub mean
        0.01, -30 ## sub sd
    )
    for(i in 1:p) {
        ui[2 * (i - 1) + 1, i] <- 1
        ui[2 * (i - 1) + 2, i] <- -1
    }
    ##
    if (verbose) {
        ll_f(theta, pheno = pheno, verbose = TRUE, age_size = age_size)
    }
    results <- constrOptim(
        theta = theta, f = ll_f, grad = NULL, ui = ui, ci = ci, mult = -1, pheno = pheno, verbose = verbose, age_size = age_size
    )
    if (calculate_ci) {
        best_theta <- results$par
        suppressWarnings(rm(ci_results_all))
        ci_results_all <- lapply(c(TRUE, FALSE), function(optimize_below) {
            mclapply(1:length(best_theta), mc.cores = nCores, function(i_param) {
                ci_results <- calculate_ci(
                    best_theta = best_theta,
                    i_param = i_param,
                    pheno = pheno,
                    age_size = age_size,
                    ci = ci,
                    ui = ui,
                    n_ci_iterations = n_ci_iterations,
                    optimize_below = optimize_below,
                    verbose = ci_verbose
                )
                return(ci_results)
            })
        })
        ## build single CI matrix here
        print(ci_results_all)
        ci_matrix <- cbind(
            lower = sapply(ci_results_all[[1]], function(x) x[["ci_value"]]),
            upper = sapply(ci_results_all[[2]], function(x) x[["ci_value"]])
        )
    } else {
        ci_matrix <- NULL
        ci_results_all <- NULL
    }
    return(
        list(
            theta = theta,
            results = results,
            ci_matrix = ci_matrix,
            ci_results_all = ci_results_all    
        )
    )
}

optimize_against_fixed_value <- function(fixed_value, fixed_name, i_param, best_theta, ui, ci, age_size, pheno) {
    theta <- best_theta[-match(fixed_name, names(best_theta))]
    ui_local <- ui[-c(2 * i_param - 1:0), -i_param]
    ci_local <- ci[-c(2 * i_param - 1:0)]
    mult <- -1
    ## werwer
    local_results <- constrOptim(theta = theta, f = ll_f, grad = NULL, ui = ui_local, ci = ci_local, mult = mult, pheno = pheno, verbose = FALSE, age_size = age_size, fixed_name = fixed_name, fixed_value = fixed_value)
    ## local_results
    ## ll_f(
    ## theta,
    ## pheno,
    ## mult = mult,
    ## verbose = TRUE,
    ## age_size = age_size,
    ## fixed_name = fixed_name,
    ## fixed_value = fixed_value
    ## )
    par <- local_results$par
    par[fixed_name] <- fixed_value
    return(list(ll = mult * local_results$value, par = par, convergence = local_results$convergence))
}

## where optimize_below = TRUE  means do CI from lower bound to best value
## where optimize_below = FALSE means do CI from upper bound to best value
calculate_ci <- function(best_theta, i_param, pheno, age_size, ci, ui, n_ci_iterations = 10, optimize_below = TRUE, verbose = FALSE) {
    fixed_name <- names(best_theta)[i_param]
    fixed_value <- best_theta[i_param]
    if (verbose) {
        print("Start CI")
        print(paste0("optimize_below = ", optimize_below))
        print(paste0("param = ", fixed_name))
    }
    ## 
    ll <- ll_f(best_theta, pheno = pheno, verbose = FALSE, age_size = age_size)
    ll_store <- array(0, c(1, 4))
    index <- 1
    ll_store[index, ] <- c(best_theta[i_param], ll, 1, 1)
    ##
    ## do 
    ##
    if (optimize_below) {
        lower <- ci[2 * i_param - 1]
        upper <- fixed_value
        val <- lower
    } else {
        lower <- fixed_value
        upper <- -1 * ci[2 * i_param - 0]
        val <- upper
    }
    ##
    r <- optimize_against_fixed_value(
        fixed_value = val,
        fixed_name = fixed_name,
        i_param = i_param,
        best_theta = best_theta,
        ui = ui,
        ci = ci,
        age_size = age_size,
        pheno = pheno
    )
    x <- (2 * (ll - r[["ll"]])) < qchisq(p = 0.95, df = 1)
    ll_store <- rbind(ll_store, c(val, r[["ll"]], x, r$convergence))
    ll_store <- ll_store[order(ll_store[, 1]), ]
    ## now, if it is inside already, screwed!
    if (x == 1) {
        print("failed!")
        ci_value <- val
        names(ci_value) <- fixed_name
        print(ll_store)
        return(list(ll_store = ll_store, ci_value = ci_value))    
    }
    ##    
    for(it in 1:n_ci_iterations) {
        val <- (lower + upper) / 2
        if (verbose) {
            print(ll_store)
            print(paste0("Trying:", val, " using lower=", lower, " and upper=", upper))
            ##theta <- best_theta
            ##theta[fixed_name] <- val
            ##ll <- ll_f(theta, pheno = pheno, verbose = FALSE, age_size = age_size)
        }
        r <- optimize_against_fixed_value(
            fixed_value = val,
            fixed_name = fixed_name,
            i_param = i_param,
            best_theta = best_theta,
            ui = ui,
            ci = ci,
            age_size = age_size,
            pheno = pheno
        )
        index <- index + 1
        inside <- (r[["ll"]] - ll) > qchisq(p = 0.95, df = 1)
        ll_store <- rbind(ll_store, c(val, r[["ll"]], (2 * (ll - r[["ll"]])) < qchisq(p = 0.95, df = 1), r$convergence))
        ## now, choose lowest 0, and highest 1        
        ll_store <- ll_store[order(ll_store[, 1]), ]
        if (optimize_below) {
            lower <- ll_store[which.max(ll_store[, 3] == 1) - 1, 1]
            upper <- ll_store[which.max(ll_store[, 3] == 1), 1]
        } else {
            lower <- ll_store[which.max(ll_store[, 3] == 0) - 1, 1]
            upper <- ll_store[which.max(ll_store[, 3] == 0), 1]
        }
        ## now, re-evaluate given boundary points
    }
    ##
    if (optimize_below) {
        ci_value <- ll_store[which.max(ll_store[, 3] == 1) - 1, 1] ## first one outside (is a 0)
    } else {
        ci_value <- ll_store[which.max(ll_store[, 3] == 0), 1] ## first one outside (is a 0)
    }
    if (verbose) {
        print("Done CI")
    }
    return(list(ll_store = ll_store, ci_value = ci_value))
}
    



try_optimizing <- function() {

    ## can I fit a full model?
    ## if independent,
    ## K_scz = 20%. then "K_sub" above should be (1 / (1 - K_scz))? 

    ## try simulating multiple times
    ## perpetually underestimate
    library("parallel")
    out <- mclapply(1:4, mc.cores = nCores, function(i_repeat) {
        set.seed(i_repeat)        
        pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
        model_out <- fit_model(pheno, verbose = FALSE, n_ci_iterations = 5, calculate_ci = TRUE, nCores = 10, ci_verbose = FALSE)
        save(pheno, model_out, file = file.path(results_dir, paste0("learn.", i_repeat, ".RData")))
        return(model_out)
    })

    lapply(out, function(x) x[["ci_matrix"]])

    ## can then do things like check CI on each parameter!
    ## do above but do it ~100 times?
    ## do CI estimate?
    ## for each, calculate 

    out

    ## add expected values
    cbind(
        default_params,
        round(out, 3)
    )
    ## if no correlation, then MLE approach correct

    ## theta <- out$theta    
    ## ll_f(theta = theta, pheno = pheno)
    ## theta["age_shape1"] <- 20
    ## ll_f(theta = theta, pheno = pheno)    
    
    ## have phenoS2
    ## sanitize a smidgen
}


fit_model_real_data <- function(pheno_file) {
    source(file.path(R_dir, "functions.R"))
    source(file.path(R_dir, "analysis_functions.R"))
    pheno <- read.csv(pheno_file)
    overlap_samples <- as.character(read.table(clozuk_overlap)[, 1])
    pheno <- pheno[-match(overlap_samples, pheno[, "IID"]), ]
    ## 
    pheno$binary_VIQ_decline <- as.integer(pheno[, "VIQ_deltazFL"] > 0.5)
    pheno$binary_FSIQ_decline <- as.integer(pheno[, "FSIQ_deltazFL"] > 0.5)
    ## I think this is just for certain histograms
    phenoS2 <- make_pheno_with_ordered_group2018(pheno)
    pheno2 <- data.frame(
        age = phenoS2[, "maxassessmentage"],
        group2018 = phenoS2[, "group2018"],
        PRS_scz = phenoS2[, "PRS_PGC_2014_SCZ"], 
        stringsAsFactors = FALSE
    )
    to_remove <-
        is.na(pheno2[, "age"]) | 
        (pheno2[, "group2018"] == "PutativeSubthreshold") & (pheno2[, "age"] >= 25) |
        (pheno2[, "group2018"] == "PutativeControl") & (pheno2[, "age"] >= 25) |
        (pheno2[, "group2018"] == "Control") & (pheno2[, "age"] < 25)
    pheno <- pheno2[!to_remove, ]
    real_pheno <- pheno
    ##
    out <- fit_model(pheno = real_pheno, calculate_ci = TRUE, nCores = 4, n_ci_iterations = 15, ci_verbose = TRUE)
    results <- cbind(out$results$par,     out$ci_matrix)
    colnames(results)[1] <- "estimate"
    results <- round(results, 2)
    write.table(
        results,
        file = file.path(results_dir, "parameter.estimates.csv"),
        row.names = TRUE,
        col.names = TRUE,
        sep = ",",
        quote = FALSE
    )
    ## 
    estimated_model_params <- get_default_model_params()
    estimated_model_params$age_size <- 70
    estimated_model_params$N <- nrow(real_pheno)
    for(name in rownames(results)) {
        estimated_model_params[[name]] <- results[name, "estimate"]
    }
    ## age fit, pretty good
    png(file = file.path(results_dir, "validate.agedist.png"), height = 6, width = 6, units = "in", res = 100)
    par(mfrow = c(2, 1))
    plot(sort(pheno[, "age"]), ylim = c(5, 80), main = "Real distribution", ylab = "Age")
    x <- rbetabinom.ab(
        n = nrow(real_pheno),
        size = estimated_model_params$age_size,
        shape1 = estimated_model_params$age_shape1,
        shape2 = estimated_model_params$age_shape2
    ) + 5
    ## 
    plot(sort(x), ylim = c(5, 80), main = "Example simulated distribution", ylab = "Age")
    dev.off()
    ##plot(log10(dbeta(x = seq(0, 1, length.out = 100), shape1 = 0.001, shape2 = 20)))
    ##plot(sort(pheno[, "age"]))
    ## plot coloured breakdown per-age, expectation from fit, vs reality
    ## do in 5 year bins?
    cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    get_counts <- function(pheno, ylab) {
        b <- ceiling(pheno[, "age"] / 5)        
        counts <- sapply(1:20, function(i) {
            x <- pheno[b == i, "group2018"]
            c(sum(x == "Case_SSD"), sum(x == "PutativeSubthreshold"), sum(x == "PutativeControl"), sum(x == "Control"))
        })
        ## plot
        plot(x = 0, y = 0, xlim = c(1, 20), ylim = c(0, 1), axes = FALSE, col = "white", xlab = "", ylab = ylab)
        for(i in 1:20) {
            a <- counts[, i]
            y <- cumsum(a) / sum(a)            
            for(j in 1:4) {
                rect(xleft = i - 0.5, xright = i + 0.5, ybottom = c(0, y), c(y, 1), col = cbPalette[2:5])
            }
        }
        par(xpd = TRUE)        
        graphics::text(x = 1:20, y = 1.2, labels = paste0("N = ", colSums(counts)), srt = 90)
        age_bins <- paste0("[", paste0(0:19 * 5 + 1), "-", paste0(0:19 * 5 + 5), "]")
        graphics::text(x = 1:20, y = -0.2, labels = age_bins, srt = 90)
    }
    png(file = file.path(results_dir, "validate.groups.png"), height = 8, width = 8, units = "in", res = 100)    
    par(mfrow = c(2, 1))
    get_counts(pheno, ylab = "Real data")
    ## now plot 
    sim_pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
    get_counts(sim_pheno, "Simulated data")
    dev.off()
    return(list(results = results, estimated_model_params = estimated_model_params))
}


get_validation_model_params <- function() {
    ##
    model_params <- get_default_model_params()
    model_params$N <- 1500 ## choose larger N for testing
    model_params$age_size <- 70
    model_params$K_scz <- 0.2
    model_params$K_sub <- 0.2
    r_e <- model_params$r_e
    r_e[2, 1] <- r_e[1, 2] <- 0
    model_params$r_e <- r_e
    r_g <- model_params$r_g
    r_g[2, 1] <- r_g[1, 2] <- 0
    model_params$r_g <- r_g
    ## here!    
    model_params$K_scz <- 0.25
    model_params$K_sub <- 0.30
    model_params$age_shape1 <- 3
    model_params$age_shape2 <- 5
    model_params$scz_mean_age <- 20
    model_params$scz_sd_age <- 2
    model_params$sub_mean_age <- 10
    model_params$sub_sd_age <- 2
    ## 
    ## default_params <- c(
    ##     model_params$K_scz,
    ##     model_params$K_sub,
    ##     model_params$age_shape1,
    ##     model_params$age_shape2,
    ##     model_params$scz_mean_age,
    ##     model_params$scz_sd_age,
    ##     model_params$sub_mean_age,
    ##     model_params$sub_sd_age
    ## )
    return(model_params)
}
    


get_default_model_params <- function() {
    ## dmp = default model params
    dmp <- list(
        N = 1500,
        K_scz = 0.2,         ## lifetime prevalence, ~ informed by data
        K_sub = 0.2,          ## lifetime prevalence, ~ informed by data
        K_rb = 0.2,           ## random binary phenotype, to inform sub-threshold
        h2_g = c(        ## everything split into known, unknown
            scz.known = 0.08, ## from our data
            scz.unknown = 0.46 - 0.08, ## gazal 2017 (Alkes price)
            sub.known = 0.08, ## not used
            sub.unknown = 0.46 - 0.08, ## assumed
            iq.known = 0.04, ## from davies 2018, prediction
            iq.unknown = 0.25 - 0.04, ## from davies 2018, table 1
            iqDecline.known = 0.20, ## irrelevant
            iqDecline.unknown = 0.40 - 0.20, ## unknown
            rb.known = 0.1, ## random binary phenotype, no genetic correlation
            rb.unknown = 0.2 ## random binary phenotype, no genetic correlation
        )
    )
    ## dmp <- append(
    ##     dmp,
    ##     list(
    ##         scz_z_thresh = qnorm(p = 1 - dmp$K_scz), ## binary thresholds
    ##         sub_z_thresh = qnorm(p = 1 - dmp$K_sub), 
    ##         rb_z_thresh = qnorm(p = 1 - dmp$K_sub)
    ##     )
    ## )
    short_pheno_names <- c("scz", "sub", "iq", "iqDecline", "rb")
    n_pheno <- length(short_pheno_names)
    ##
    dmp <- append(
        dmp,
        list(
            n_pheno = n_pheno
        )
    )
    ##
    r_g <- array(0, c(n_pheno, n_pheno))
    diag(r_g) <- 1 ## by definition
    i <- 1; j <- 2; r_g[i, j] <- r_g[j, i] <-  0.80  ## scz, sub
    i <- 1; j <- 3; r_g[i, j] <- r_g[j, i] <- -0.234 ## scz, iq davies 2018 supp data 9
    i <- 1; j <- 4; r_g[i, j] <- r_g[j, i] <-  0.00  ## scz, iqDecline
    i <- 2; j <- 3; r_g[i, j] <- r_g[j, i] <-  0.00  ## sub, iq
    i <- 2; j <- 4; r_g[i, j] <- r_g[j, i] <-  0.00  ## sub, iqDecline
    i <- 3; j <- 4; r_g[i, j] <- r_g[j, i] <-  0.80  ## iq, iqDecline
    colnames(r_g) <- short_pheno_names
    rownames(r_g) <- colnames(r_g)
    dmp <- append(dmp, list(r_g = r_g))
    ## correlation in environmental noise
    r_e <- array(NA, c(n_pheno, n_pheno))
    r_e[] <- 0.5 ## set rest to 0.5
    diag(r_e) <- 1 ## by definition
    ## i <- 1; j <- 2; r_e[i, j] <- r_e[j, i] <- 0.5 ## scz, sub
    ## i <- 1; j <- 3; r_e[i, j] <- r_e[j, i] <- 0.5 ## scz, iq
    ## i <- 1; j <- 4; r_e[i, j] <- r_e[j, i] <- 0.5 ## scz, iqDecline
    ## i <- 2; j <- 3; r_e[i, j] <- r_e[j, i] <- 0.5 ## sub, iq
    ## i <- 2; j <- 4; r_e[i, j] <- r_e[j, i] <- 0.5 ## sub, iqDecline
    ## i <- 3; j <- 4; r_e[i, j] <- r_e[j, i] <- 0.5 ## scz, iqDecline
    colnames(r_e) <- short_pheno_names
    rownames(r_e) <- colnames(r_e)
    ##
    dmp <- append(dmp, list(r_e = r_e))
    ##age_min = 10,
    ##age_max = 70,
    dmp <- append(
        dmp,
        list(
            age_shape1 = 1.5,
            age_shape2 = 3.5,
            scz_mean_age = 20,
            scz_sd_age = 10,
            sub_mean_age = 10,
            sub_sd_age = 2,
            rb_mean_age = 16,
            rb_sd_age = 2,
            age_size = 80
        )
    )
    return(dmp)
}




## check that variables are sufficiently similar
check1 <- function(a, b, tol = 0.05) {
    if (abs(a - b) > tol) {
        print(a)
        print(b)
        stop("bad math")
    }
}
## check things are not both true
check2 <- function(a, b) {
    if (sum(a & b) > 0) {
            stop("bad code")
    }
}


if (1 == 0) {

    o = array(0, c(5, 5))
    for(i in 1:5) {
        for(j in 1:5) {
            o[i, j] = cor(E[, i], E[, j])
        }
    }

}




simulate_full <- function(model_params = get_default_model_params(), do_checks = TRUE) {

    ## inject into local environment
    for (i in seq_along(model_params)) {
        assign(names(model_params)[i], model_params[[i]])
    }
    
    ## determine thresholds
    scz_z_thresh <- qnorm(p = 1 - K_scz)
    sub_z_thresh <- qnorm(p = 1 - K_sub)
    rb_z_thresh <- qnorm(p = 1 - K_sub)

    ## simulate according to r_g
    SigmaG <- array(NA, c(n_pheno, n_pheno))
    SigmaE <- array(NA, c(n_pheno, n_pheno))
    for(i in 1:n_pheno) {
        for(j in 1:n_pheno) {
            n1 <- colnames(r_g)[i]
            n2 <- colnames(r_g)[j]
            h2_g1 <- sum(h2_g[grep(paste0(n1, "."), names(h2_g), fixed = TRUE)])
            h2_g2 <- sum(h2_g[grep(paste0(n2, "."), names(h2_g), fixed = TRUE)])
            SigmaG[i, j] <- r_g[i, j] * sqrt(h2_g1) * sqrt(h2_g2)
            SigmaE[i, j] <- r_e[i, j] * sqrt(1 - h2_g1) * sqrt(1 - h2_g2)
        }
    }

    G <- MASS::mvrnorm(
        n = N,
        mu = rep(0, ncol(r_e)),
        Sigma = SigmaG
    )
    E <- MASS::mvrnorm(
        n = N,
        mu = rep(0, ncol(r_e)),
        Sigma = SigmaE
    )
    colnames(G) <- colnames(r_g)
    colnames(E) <- colnames(r_e)    

    for(i in 1:n_pheno) {
        G_3 <- G[, i]
        p <- colnames(G)[i]
        s1 <- h2_g[paste0(p, ".known")]
        s2 <- h2_g[paste0(p, ".unknown")]
        mean <- G_3 * s1 / (s1 + s2)
        var <- (s1 * s2) / (s1 + s2)
        G_1 <- rnorm(N, mean = mean, sd = sqrt(var))
        G_2 <- G_3 - G_1
        new_cols <- paste0(p, c(".known", ".unknown"))
        if (sum(new_cols %in% colnames(G)) == 0) {
            G <- cbind(G, G_1, G_2)
            colnames(G)[ncol(G) + -1:0] <- new_cols
        }
    }

    ## do checks here
    if (do_checks) {
        ## check everything
        for(i in 1:ncol(E)) {
            n1 <- colnames(r_e)[i]
            h2_g1 <- h2_g[grep(paste0(n1, "."), names(h2_g), fixed = TRUE)]
            check1(var(G[, i]), sum(h2_g1))
            check1(var(E[, i]), 1 - sum(h2_g1))
        }
        ## correlations
        for(i in 1:4) {
            for(j in 1:4) {
                check1(cor(G[, i], G[, j]), r_g[i, j])                
                check1(cor(E[, i], E[, j]), r_e[i, j])
            }
        }
        for(i in 1:length(h2_g)) {
            check1(var(G[, names(h2_g)[i]]), h2_g[i])
        }
    }

    ## OK, build!
    Y_scz <- G[, "scz"] + E[, "scz"]
    Y_sub <- G[, "sub"] + E[, "sub"]
    Y_rb <- G[, "rb"] + E[, "rb"]
    Y_iq <- G[, "iq"] + E[, "iq"]
    Y_iqDecline <- G[, "iqDecline"] + E[, "iqDecline"]
    if (do_checks) {
        check1(var(Y_scz), 1)
        check1(var(Y_sub), 1)
        check1(var(Y_rb), 1)
        check1(var(Y_iq), 1)
        check1(var(Y_iqDecline), 1)        
    }

    ## note - 25 means control. 24 and under means putative
    
    ## build binary phenotypes
    T_scz_age <- rnorm(N, mean = scz_mean_age, sd = scz_sd_age) ## threshold age, if develop SCZ
    T_sub_age <- rnorm(N, mean = sub_mean_age, sd = sub_sd_age) ## sub-threshold age
    ## meh
    raw_age <- rbetabinom.ab(n = N, size = age_size, shape1 = age_shape1, shape2 = age_shape2)
    ##print(range(raw_age))
    age <- raw_age + 5
    ##print(range(age))
    ##print(mean(age))
    ##print(age_size)
    ## age <- age_min + (age_max - age_min) * rbeta(N, shape1 = age_shape1, shape2 = age_shape2)

    ## shape2 higher pushes leftwards
    Z_scz <- (scz_z_thresh < Y_scz) & (T_scz_age < age) ## has SCZ, above age
    Z_sub <- (Z_scz == 0) & (sub_z_thresh < Y_sub) & (T_sub_age < age) & (age < 25) ## no SCZ, has sub-threshold
    Z_putcont <- (Z_scz == 0) & (Z_sub == 0) & (age < 25)    
    Z_cont <- (Z_scz == 0) & (Z_sub == 0) & (25 <= age)
    ## random binary phenotype
    T_rb_age <- rnorm(N, mean = rb_mean_age, sd = rb_sd_age) ## random binary age
    Z_rb <- (Z_scz == 0) & (age < 25) & (rb_z_thresh < Y_rb) & (T_rb_age < age) ## no scz, has random-binary phenotype
    Z_putcont_rb <- (Z_scz == 0) & (Z_rb == 0) & (age < 25)
    Z_cont_rb <- (Z_scz == 0) & (Z_rb == 0) & (25 <= age)
    ## does it ever meet phenotype?
    ever_scz <- (scz_z_thresh < Y_scz) 
    ever_sub <- (sub_z_thresh < Y_sub) 
    ever_rb <- (rb_z_thresh < Y_rb) 

    
    if (do_checks) {
        ## normal
        check2(Z_scz, Z_sub)
        check2(Z_scz, Z_putcont)
        check2(Z_scz, Z_cont)
        check2(Z_sub, Z_putcont)
        check2(Z_sub, Z_cont)
        check2(Z_putcont, Z_cont)
        ## random binary
        check2(Z_scz, Z_rb)
        check2(Z_scz, Z_putcont_rb)
        check2(Z_scz, Z_cont_rb)
        check2(Z_rb, Z_putcont_rb)
        check2(Z_rb, Z_cont_rb)
        check2(Z_putcont_rb, Z_cont_rb)
    }

    group2018 <- array("", N)
    group2018[Z_scz] <- "Case_SSD"
    group2018[Z_sub] <- "PutativeSubthreshold"
    group2018[Z_putcont] <- "PutativeControl"
    group2018[Z_cont] <- "Control"

    ## do with random binary
    ## KEEP SAME NAMES for simplicity
    group2018_rb <- array("", N)
    group2018_rb[Z_scz] <- "Case_SSD"
    group2018_rb[Z_rb] <- "PutativeSubthreshold"
    group2018_rb[Z_putcont_rb] <- "PutativeControl"
    group2018_rb[Z_cont_rb] <- "Control"

    ## so far, keep un-related
    sex <- sample(c(1, 2), N, replace = TRUE)

    ## define binaryIQ decline
    iqDeclineBinary <- as.integer((Y_iqDecline) > 0.5)

    ## TRANSFORM NORMAL -> QUESTION
    ## make new alternative subthresold quantitative phenotype through approximate quantiles
    subQ2 <- round(qexp(pnorm(Y_sub), rate = 0.2238))
    
    pheno <- data.frame(
        group2018 = group2018,
        group2018_rb= group2018_rb,
        PRS_scz = G[, "scz.known"],
        PRS_iq = G[, "iq.known"],
        age = age,
        sex = sex,
        iq = Y_iq,
        sub = Y_sub,
        subQ = round(5 ** Y_sub),
        subQ2 = subQ2,
        iqDecline = Y_iqDecline,
        full_scz = G[, "scz"],
        ever_scz = ever_scz,
        ever_sub = ever_sub,
        ever_rb = ever_rb        
    )

    return(
        list(
            pheno = pheno
        )
    )
    
}    





if (1 == 0) {


    c <- "group2018"
    c <- "group2018_rb"
    gg <- c("Case_SSD", "Control", "PutativeControl", "PutativeSubthreshold")
    sapply(gg, function(g) {
        w <- pheno[, c] == g
        return(mean(pheno[w, "PRS_scz"]))
    })
    mean(pheno[pheno[, c] == "PutativeControl" | pheno[, c] == "Control", "PRS_scz"])
    ## so very different between PutativeSubthreshold and controls? when effect exists
    ## correct comparison indeed seems to be combining them

    ## can I build components, then sum?
    nPheno <- length(h2_g)
    SigmaG <- array(NA, c(nPheno, nPheno))
    diag(SigmaG) <- h2_g
    pheno_names <- unique(sapply(strsplit(names(h2_g), ".", fixed = TRUE), I)[1, ])
    rownames(SigmaG) <- names(h2_g)
    colnames(SigmaG) <- rownames(SigmaG)


    ## off-target
    for(i1 in 1:nPheno) {
        for(i2 in 1:nPheno) {
            v1 <- h2_g[i1]
            v2 <- h2_g[i2]
            n1 <- strsplit(names(v1), ".", fixed = TRUE)[[1]][1]
            o1 <- strsplit(names(v1), ".", fixed = TRUE)[[1]][2]            
            n2 <- strsplit(names(v2), ".", fixed = TRUE)[[1]][1]
            o2 <- strsplit(names(v2), ".", fixed = TRUE)[[1]][2]            
            sigma_pheno1 <- sqrt(sum(h2_g[grep(n1, names(h2_g))]))
            sigma_pheno2 <- sqrt(sum(h2_g[grep(n1, names(h2_g))]))
            ## if known / unknown of same phenotype, set to 0 by default
            ##if ((o1 != o2) & (n1 == n2)) {
            ##    SigmaG[names(v1), names(v2)] <- 0
            ##} else {
            SigmaG[names(v1), names(v2)] <- sqrt(v1) * sqrt(v2) * r_g[n1, n2]
            ##}
        }
    }
    
    ##
    G <- MASS::mvrnorm(
        n = N,
        mu = rep(0, 4),
        Sigma = SigmaG[1:4, 1:4]
        )
    
    colnames(G) <- names(h2_g)

    if (do_checks) {
        ## check everything
        for(i in 1:ncol(G)) {
            check1(var(G[, i]), h2_g[i])
        }
        ## constructed variables
        for(i in 1:4) {
            check1(
                var(G[, 2 * i - 1] + G[, 2 * i]),
                sum(h2_g[2 * i -1:0])
            )
            for(j in 1:4) {
                check1(
                    cor(
                        G[, 2 * i - 1] + G[, 2 * i],
                        G[, 2 * j - 1] + G[, 2 * j]
                    ),
                    r_g[i, j]
                )
            }
        }
    }

}



simulate_group2018 <- function(do_checks = TRUE) {

    ## see https://en.wikipedia.org/wiki/Covariance#Properties
    ## and https://en.wikipedia.org/wiki/Covariance_and_correlation

    ## latent factors
    G_l1 <- rnorm(N, mean = 0, sd = 1) ## known SCZ
    G_l2 <- rnorm(N, mean = 0, sd = 1) ## unknown SCZ
    G_l3 <- rnorm(N, mean = 0, sd = 1) ## specific to sub-threshold
    ## genetic factors
    G_scz <- sqrt(h2_g_scz_known) * G_l1 + sqrt(h2_g_scz_unknown) * G_l2
    ## let X be (G_scz / h2_g_scz) ~ N(0, 1)
    ## then we are interested in variables c and d such that 
    ## then cor(aX, cX + dV) * h_scz * h_sub = cov(aX, cX + dV)
    ## then r * h_sub * h_scz = c * h_scz
    ## then c = h_sub * r
    ## then h2_sub = var(cX + dV) = c^2 + d ^ 2
    ## then d = sqrt(h2_sub - c^2) = sqrt(h2_sub (1 - r2))
    ## make sure to use X defined properly in above and it should work, subject to constraints between correlation and heritabilities
    ## liabilities
    G_sub <-
        sqrt(h2_g_sub) * sqrt(r2_g_scz_sub) * (G_scz / sqrt(h2_g_scz)) +
        sqrt(h2_g_sub * (1 - r2_g_scz_sub)) * G_l3
    ## check math checks out
    if (do_checks) {
        check1(var(G_sub), h2_g_sub)
        check1(cor(G_scz, G_sub) ** 2, r2_g_scz_sub)
        check1(var(G_scz), h2_g_scz)
    }

    ## environmental factors
    eps1 <- rnorm(N, mean = 0, sd = 1) ## joint 
    eps2 <- rnorm(N, mean = 0, sd = 1) ## scz specific
    eps3 <- rnorm(N, mean = 0, sd = 1) ## sub specific
    ## want environmental correlation of r2_e_scz_sub
    ## here X, Y and Z are all N(0, 1)
    ## cor(aX + bY, cX + dZ) * h_e_scz * h_e_sub = cov(aX + bY, cX + dZ) = a * c
    ## so r_e * h_e_scz * h_e_sub = a * c
    ## a^2 + b^2 = h2_e_scz
    ## c^2 + d^2 = h2_e_sub
    ## let  a = r_e * h_e_scz
    ##      c = h_e_sub
    ## then b = sqrt(h2_e_scz * (1 - r_e * h_e_scz))
    ##      d = 0
    eps_scz <- sqrt(h2_e_scz) * (
        sqrt(r2_e_scz_sub) * eps1 +
        sqrt(1 - r2_e_scz_sub) * eps2
    )
    eps_sub <- sqrt(h2_e_sub) * (
        eps1
    )
    ## environmental terms
    if (do_checks) {
        check1(cor(eps_scz, eps_sub) ** 2, r2_e_scz_sub)
        check1(var(eps_scz), h2_e_scz)
        check1(var(eps_sub), h2_e_sub)
    }
    
    
    ## now, build rest
    Y_scz <- G_scz + eps_scz
    Y_sub <- G_sub + eps_sub
    if (do_checks) {
        check1(var(Y_scz), 1)
        check1(var(Y_sub), 1)
    }
    T_scz_age <- rnorm(N, mean = scz_mean_age, sd = scz_sd_age) ## threshold age, if develop SCZ
    T_sub_age <- rnorm(N, mean = 16, sd = 2) ## sub-threshold age
    age <- age_min + (age_max - age_min) * rbeta(N, shape1 = age_shape1, shape2 = age_shape2) ## shape2 higher pushes leftwards
    
    Z_scz <- (scz_z_thresh < Y_scz) & (T_scz_age < age) ## has SCZ, above age
    Z_sub <- (sub_z_thresh < Y_sub) & (T_sub_age < age) & (Z_scz == 0) ## no SCZ, has sub-threshold
    Z_putcont <- (Z_scz == 0) & (Z_sub == 0) & (age < 25)
    Z_cont <- (Z_scz == 0) & (Z_sub == 0) & (25 < age)

    if (do_checks) {
        check2(Z_scz, Z_sub)
        check2(Z_scz, Z_putcont)
        check2(Z_scz, Z_cont)
        check2(Z_sub, Z_putcont)
        check2(Z_sub, Z_cont)
        check2(Z_putcont, Z_cont)
    }

    group2018 <- array("", N)
    group2018[Z_scz] <- "Case_SSD"
    group2018[Z_sub] <- "PutativeSubthreshold"
    group2018[Z_putcont] <- "PutativeControl"
    group2018[Z_cont] <- "Control"

    ## so far, keep un-related
    sex <- sample(c(1, 2), N, replace = TRUE)
    
    pheno <- data.frame(
        group2018 = group2018,
        PRS_scz = G_l1,
        age = age,
        sex = sex
    )

    return(
        list(
            pheno = pheno
        )
    )
    
}



## head to head statistics
analyze_pheno_for_Aim2A <- function(pheno, age = "age", sex = "sex", group2018_name = "group2018", which = NULL, prs_colname = "PRS_scz") {
    p_PRS <- NULL
    if (is.null(which)) {
        which <- 1:7
    } 
    for(i in which) {
        var1 <- c("Case_SSD", "Control",         "PutativeSubthreshold", "Case_SSD",             "PutativeSubthreshold", "PutativeSubthreshold", "Case_SSD")[i]
        ## next two together
        var2 <- c("Control" , "PutativeControl", "Control"             , "PutativeSubthreshold", "PutativeControl",      "Control",              "Control")[i]
        var3 <- c("",         "",                "PutativeControl"     , "",                     "",                     "",                     "PutativeControl")[i]
        bin <- NA
        bin[pheno[, group2018_name] == var1] <- 1
        bin[pheno[, group2018_name] == var2] <- 0
        if (var3 != "") {
            bin[pheno[, group2018_name] == var3] <- 0
            var2 <- paste0(var2, "&", var3)
        }
        ##
        g <- glm(bin ~ pheno[, prs_colname], family = "binomial")        
        ##
        if (!g$converged) {
            print(i)
            print("not converged")
            g <- glm(bin ~ pheno[, prs_colname], family = "binomial")
        }
        c1 <- coefficients(summary(g))
        ##   
        p_PRS <- c(p_PRS, c1[2, 4])
        names(p_PRS)[length(p_PRS)] <- paste0(var1, "_vs_", var2)
        ##print("-------------")
        ##print(p_PRS[length(p_PRS)])
        ##print(sum(!is.na(bin)))
    }
    return(p_PRS)
}


## head to head statistics
analyze_pheno_for_Aim2AB <- function(
    pheno,
    subpheno = "subQ",
    age = "age",
    sex = "sex",
    group2018_name = "group2018",
    include_aim2A = FALSE
) {
    p_PRS <- NULL
    ## from simulations, should substantially outnumber...
    ## just calculate quantitative subthreshold among those with subthreshold / putative control?
    c1 <- NA
    c2 <- NA
    for(i in 1:3) {
        var1 <- c("PutativeSubthreshold", "PutativeControl", "PutativeSubthreshold")[i]
        var2 <- c(                    "",                "",      "PutativeControl")[i]
        bin <- array(0, nrow(pheno))
        bin[pheno[, group2018_name] == var1] <- 1
        bin[pheno[, group2018_name] == var2] <- 1
        phenoL <- pheno[bin == 1, ]
        g1 <- lm(as.formula(paste0(subpheno, " ~ PRS_scz")),  data = phenoL)
        c1 <- coefficients(summary(g1))
        check1 <- tryCatch({
            g2A <- lm(as.formula(paste0(subpheno, " ~ PRS_scz + ", group2018_name)),  data = phenoL)            
            0
        },
        error = function(e) 1
        )
        check2 <- tryCatch({
            g2B <- lm(as.formula(paste0(subpheno, " ~ PRS_scz")),  data = phenoL)
            0
        },
        error = function(e) 1
        )
        if ((check1 == 0) & (check2 == 0)) {        
            c2 <- coefficients(summary(g2A))
            c3 <- coefficients(summary(g2B))            
        }
        f <- function(var1, var2, m) {
            if (var2 == "") {
                return(var1)
            } else {
                return(paste0(var1, "_with_", var2, ".", m))
            }
        }
        if (i <= 2) {
            p_PRS <- c(p_PRS, c1[2, 4])
            name <- f(var1, var2, "adjust")
            names(p_PRS)[length(p_PRS)] <- name
            ## N
            p_PRS <- c(p_PRS, nrow(phenoL))
            names(p_PRS)[length(p_PRS)] <- paste0("N_", name)
        } else {
            p_PRS <- c(p_PRS, c2[2, 4])
            name <- f(var1, var2, "adjust")
            names(p_PRS)[length(p_PRS)] <- name
            ## N
            p_PRS <- c(p_PRS, nrow(phenoL))
            names(p_PRS)[length(p_PRS)] <- paste0("N_", name)
            ##            
            p_PRS <- c(p_PRS, c3[2, 4])
            name <- f(var1, var2, "noadjust")
            names(p_PRS)[length(p_PRS)] <- name
            ## N
            p_PRS <- c(p_PRS, nrow(phenoL))
            names(p_PRS)[length(p_PRS)] <- paste0("N_", name)
        }
    }
    ##
    ## also, get p-value for existence of difference between groups
    ##
    if (include_aim2A) {
        p_PRS <- c(
            p_PRS,
            analyze_pheno_for_Aim2A(pheno, which = 3)["PutativeSubthreshold_vs_Control&PutativeControl"]
        )
    }
    ##
    ##
    ##
    return(p_PRS)
}


## should I include age here?
analyze_pheno_for_Aim2B <- function(pheno, age = "age", sex = "sex") {
    ## do tests here
    s1 <- summary(lm(iq ~ PRS_scz, data = pheno))
    s2 <- summary(lm(iq ~ PRS_iq, data = pheno))
    ## binary test
    pheno$iqDeclineBinary <- as.integer((pheno$iqDecline) > 0.5)    
    s3 <- summary(glm(iqDeclineBinary ~ PRS_scz, data = pheno, family = "binomial"))
    s4 <- summary(glm(iqDeclineBinary ~ PRS_iq, data = pheno, family = "binomial"))
    s5 <- summary(lm(iqDecline ~ PRS_scz, data = pheno))
    s6 <- summary(lm(iqDecline ~ PRS_iq, data = pheno))
    ## get p-values
    to_out <- c(
        p_iq_vs_PRS_scz =        coefficients(s1)[2, 4],
        p_iq_vs_PRS_iq =         coefficients(s2)[2, 4],
        p_iqDeclineBinary_vs_PRS_scz = coefficients(s3)[2, 4],
        p_iqDeclineBinary_vs_PRS_iq  = coefficients(s4)[2, 4],
        p_iqDecline_vs_PRS_scz = coefficients(s5)[2, 4],
        p_iqDecline_vs_PRS_iq = coefficients(s6)[2, 4]
    )
    return(to_out)
}


get_group_counts <- function(pheno, group2018_name) {
    labels <- c("Case_SSD", "Control", "PutativeControl", "PutativeSubthreshold")
    return(sapply(labels, function(label) {
        sum(pheno[, group2018_name] == label)
    }))
}

if (1 == 0) {
    
    ##
    ## this blurb looks at power
    ##
    
    file <- file.path(results_dir, "model_params.RData")    
    load(file = file)
    ## modest genetic correlation!
    model_params$r_g[, ] <- 0
    diag(model_params$r_g) <- 1
    model_params$r_g["scz", "iq"] <- model_params$r_g["iq", "scz"] <- -0.234
    model_params$N <- 200000
    rm(pheno); pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
    
    prs <- "PRS_iq"
    ## prs <- "PRS_scz"
    c(
        mean(pheno[pheno[, "group2018"] == "Case_SSD", prs])        ,
        mean(pheno[pheno[, "group2018"] == "PutativeControl", prs]),
        mean(pheno[pheno[, "group2018"] == "Control", prs])
    )
    bin1 <- array(NA, nrow(pheno))
    bin1[pheno[, "group2018"] == "Case_SSD"] <- 1
    bin1[pheno[, "group2018"] == "PutativeControl"] <- 0
    bin1[pheno[, "group2018"] == "Control"] <- 0
    ## other bin
    bin2 <- array(NA, nrow(pheno))
    bin2[pheno[, "group2018"] == "Case_SSD"] <- 1
    bin2[pheno[, "group2018"] == "Control"] <- 0
    g1a <- glm(bin1 ~ pheno[, "PRS_scz"], family = "binomial")
    g1b <- glm(bin2 ~ pheno[, "PRS_scz"], family = "binomial")    
    g2a <- glm(bin1 ~ pheno[, "PRS_iq"], family = "binomial")
    g2b <- glm(bin2 ~ pheno[, "PRS_iq"], family = "binomial")
    g3 <- lm(pheno[, "iq"] ~ pheno[, "PRS_iq"])
    ## r2 difference? iq?
    PseudoR2(g1a, "Nagelkerke") ## ~4% OK
    PseudoR2(g1b, "Nagelkerke") ## ~6% smidgen low
    PseudoR2(g2a, "Nagelkerke") ## 0.36%?? seems low?
    PseudoR2(g2b, "Nagelkerke") ## 0.4%?
    summary(g3)$r.squared ## ~4%, about right

    ## is h2 correct?
    file <- file.path(results_dir, "model_params.RData")    
    load(file = file)
    calculate_power(model_params, n_reps = 100, alpha = 0.05, group2018_name = "group2018", prs_colname = "PRS_scz")
    calculate_power(model_params, n_reps = 100, alpha = 0.05, group2018_name = "group2018", prs_colname = "PRS_iq")

}

calculate_power <- function(model_params, n_reps = 100, alpha = 0.05, group2018_name = "group2018", prs_colname = "PRS_scz") {

    results <- parallel::mclapply(1:n_reps, mc.cores = nCores, function(i) {
        check <- tryCatch({
            pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
            0
        },
        error = function(e) 1
        )
        if (check == 0) {
            return(list(
                0,
                analyze_pheno_for_Aim2A(pheno, group2018_name = group2018_name, prs_colname = prs_colname),
                get_group_counts(pheno, group2018_name)
            ))
        } else {
            return(list(1, 1))
        }
    })
    
    initialized <- FALSE
    count <- 0
    for(iRep in 1:n_reps) {
        check <- results[[iRep]][[1]]
        if (check == 0) {
            if (!initialized) {
                power <- array(0, length(results[[1]][[2]]))
                av_p <- power
                av_n <- array(0, length(results[[1]][[3]]))
                initialized <- TRUE
            }
            r <- (results[[iRep]][[2]] < alpha)
            power <- power + r
            count <- count + 1
            av_p <- av_p + results[[iRep]][[2]]
            av_n <- av_n + results[[iRep]][[3]]
        }
    }
    if (!initialized) {
        return(rep(NA, 6))
    }
    av_p <- av_p / count
    power <- power / count
    av_n <- av_n / count
    ## skeaky, get names
    pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno    
    names(power) <- names(analyze_pheno_for_Aim2A(pheno, group2018_name = group2018_name))
    names(av_p) <- names(power)
    ## 
    ## cbind(power, av_p, av_n)
    ##
    return(power)
}






power_analysis_aim2A <- function(filename, n_power_reps = 100, group2018_name = "group2018", model_params = get_default_model_params(), prs_colname = "PRS_SCZ") {

    ## does this need to be run manually?
    ## calculate "power" / max estimate for r2_g_scz_sub
    r_g <- model_params$r_g
    ## assume transitive genetic correlation
    r_g_scz_sub_range <- seq(0, 1, length.out = 21)
    power_mat <- NULL
    for(r_g_scz_sub in r_g_scz_sub_range) {
        r_g["scz", "sub"] <- r_g["sub", "scz"] <- r_g_scz_sub
        ## assume transitive
        r_g["iq", "sub"] <- r_g["sub", "iq"] <- r_g_scz_sub * r_g["scz", "iq"]
        model_params$r_g <- r_g
        message(paste0(r_g_scz_sub, ", ", date()))
        calculate_power(model_params = model_params, n_reps = n_power_reps, alpha = 0.05, group2018_name = group2018_name, prs_colname = prs_colname)
        power_mat <- rbind(
            power_mat,
            calculate_power(model_params = model_params, n_reps = n_power_reps, alpha = 0.05, group2018_name = group2018_name, prs_colname = prs_colname)
        )
        
    }
    rownames(power_mat) <- r_g_scz_sub_range
    
    ## save for now?
    save(
        power_mat,
        file = paste0(filename, ".power.stuff.RData")
    )
    
    ## colnames(power_mat) <- c(
    ##     "Case_SSD_vs_Control",
    ##     "Control_vs_PutativeControl",
    ##     "PutativeSubthreshold_vs_Control&PutativeControl",
    ##     "Case_SSD_vs_PutativeSubthreshold",
    ##     "PutativeSubthreshold_vs_PutativeControl",
    ##     "PutativeSubthreshold_vs_Control"
    ## )

    remapper <- function(x) {
        if (x == "Case_SSD_vs_Control") {
            return("CaseSSD vs Control")
        }
        if (x == "Control_vs_PutativeControl") {
            return("Control vs PutativeControl")
        }
        if (x == "PutativeSubthreshold_vs_Control&PutativeControl") {
            return("Putative Subthreshold vs Control and Putative Control")
        }
        if (x == "Case_SSD_vs_PutativeSubthreshold") {
            return("CaseSSD vs Putative Subthreshold")
        }
        if (x == "PutativeSubthreshold_vs_PutativeControl") {
            return("Putative Subthreshold vs Putative Control")
        }
        if (x == "PutativeSubthreshold_vs_Control") {
            return("Putative Subthreshold vs Control")
        }
    }
    
    ## tell
    for(suffix in c("png", "pdf")) {
        image_open(filename = filename, height = 8, width = 12, suffix = suffix)
        par(mfrow = c(2, 3))
        for(test in colnames(power_mat)) {
            plot(
                x = r_g_scz_sub_range,
                y = power_mat[, test],
                xlab = "correlation (r_g)",
                ylab = "Power",
                type = "l",
                main = gsub("_", " ", test),
                ylim = c(0, 1)
            )
        }
        dev.off()
    }

}


## here, for scz and sub
## show power as a function between scz and sub phenotype
power_plots_aim2A <- function() {

    ## on top, plot
    for(suffix in c("png", "pdf")) {
        image_open(
            filename = file.path(results_dir, paste0("aim2a.power.bothPS")),
            height = 8, width = 8, suffix = suffix
        )
        par(mfrow = c(2, 2))
        for(i_ps in 1:2) {
            prs_colname <- c("PS_sz", "PS_iq")[i_ps]
            if (i_ps == 1) {load(file.path(results_dir, "aim2a.power.power.stuff.RData")) } 
            if (i_ps == 2) {load(file.path(results_dir, "aim2a.iq.power.power.stuff.RData")) }
            ## what do I care about
            to_plot <- c(
                "Case_SSD_vs_Control&PutativeControl",
                "PutativeSubthreshold_vs_Control&PutativeControl"
            )
            ## load here!
            for(test in to_plot) {
                plot(
                    x = r_g_scz_sub_range,
                    y = power_mat[, test],
                    xlab = "correlation (r_g)",
                    ylab = "Power",
                    type = "l",
                    main = paste0(prs_colname, "\n", gsub("vs", "vs\n", gsub("_", " ", test))),
                    ylim = c(0, 1)
                )
            }
        }
        dev.off()
    }
    
    
}



calculate_aim2AB_power <- function(model_params, nRep, alpha = 0.05, include_aim2A = FALSE) {
    out <- parallel::mclapply(1:nRep, mc.cores = nCores, function(iRep) {
        check <- tryCatch({
            pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
            0
        },
        error = function(e) 1
        )
        r_sub <- NULL
        r_subQ <- NULL        
        if (check == 0) {
            ## transform to normal
            pheno$y <- qnorm(pexp(q = 0.5 + pheno[, "subQ2"], rate = 0.2238))
            ## analyze - only
            bin <- array(NA, nrow(pheno))
            bin[pheno[, group2018_name] == "PutativeSubthreshold"] <- 1
            bin[pheno[, group2018_name] == "PutativeControl"] <- 0
            phenoL <- pheno[!is.na(bin), ]
            phenoL$bin <- bin[!is.na(bin)]
            p_values <- c(
                "quantT_on_PRS_scz" = coefficients(summary(lm("y ~ PRS_scz", data = phenoL)))[2, 4],
                "quantT_on_PRS_scz_plus_bin" = coefficients(summary(lm("y ~ PRS_scz + bin", data = phenoL)))[2, 4],
                "oriBin_on_PRS_scz" = coefficients(summary(lm("bin ~ PRS_scz", data = phenoL)))[2, 4]
            )
        } else {
            p_values <- rep(NA, 3)
        }
        return(p_values)
    })
    ##
    out2 <- t(sapply(out, I))
    ##
    power_original <- colSums(out2 < alpha) / colSums(!is.na(out2))
    ## also - first two, conditional on third
    cond <- out2[, "oriBin_on_PRS_scz"] < alpha
    power_conditional <- colSums(out2[cond, , drop = FALSE] < alpha) / colSums(!is.na(out2[cond, , drop = FALSE]))
    ## 
    return(
        list(
            power_original = power_original,
            power_conditional = power_conditional
        )
    )
}


clean_aim2AB_power <- function(out, r_name = "r_sub", alpha = 0.05, include_aim2A) {
    to_out_sig <- array(0, length(out[[1]][[r_name]]))
    names(to_out_sig) <- names(out[[1]][[r_name]])
    to_out_n <- to_out_sig
    to_out_sig_cond <- to_out_sig
    to_out_n_cond <- to_out_n
    for(iRep in 1:length(out)) {
        r <- out[[iRep]][[r_name]]
        check <- out[[iRep]]$check                
        if (check == 0) {
            to_out_n <- to_out_n + 1
            to_out_sig <- to_out_sig + as.integer(r < alpha)
            if (include_aim2A) {
                if (r["PutativeSubthreshold_vs_Control&PutativeControl"] < alpha) {
                    to_out_n_cond <- to_out_n_cond + 1
                    to_out_sig_cond <- to_out_sig_cond + as.integer(r < alpha)
                }
            }
        }
    }
    to_out_power <- to_out_sig / to_out_n
    to_out_power_cond <- to_out_sig_cond / to_out_n_cond
    return(list(to_out_power = to_out_power, to_out_power_cond = to_out_power_cond))
}


## 
## power analysis for genetic correlation for versions of the quantitative phenotype
## 
power_analysis_aim2AB <- function(fileprefix, n_power_reps = 100, group2018_name = "group2018", model_params = get_default_model_params()) {

    ## does this need to be run manually?
    ## calculate "power" / max estimate for r2_g_scz_sub
    r_g_scz_sub_range <- seq(0, 1, length.out = 21)
    r_g <- model_params$r_g 
    r_g[, ] <- 0
    diag(r_g) <- 1
    power_mat_original <- NULL
    power_mat_conditional <- NULL
    for(r_g_scz_sub in r_g_scz_sub_range) {
        r_g["scz", "sub"] <- r_g["sub", "scz"] <- r_g_scz_sub
        model_params$r_g <- r_g
        message(paste0(r_g_scz_sub, ", ", date()))
        ## this is not quite what I did
        p <- calculate_aim2AB_power(model_params, n_power_reps, alpha = 0.05, include_aim2A = TRUE)
        ##
        power_mat_original <- rbind(power_mat_original, p$power_original)
        power_mat_conditional <- rbind(power_mat_conditional, p$power_conditional)        
    }

    ## save for now?
    save(
        power_mat_original,
        power_mat_conditional,
        file = paste0(fileprefix, ".power.stuff.RData")
    )

    ## tell
    for(suffix in c("png", "pdf")) {
        image_open(filename = paste0(fileprefix, ".", subpheno), height = 4, width = 8, suffix = suffix)
        par(mfrow = c(1, 2))
        tests_to_plot <- c(
            "quantT_on_PRS_scz",
            "quantT_on_PRS_scz_plus_bin"
        )
        ## 
        names_tests_to_plot <- c(
            "Quantitative SIPS vs\nPS_scz no adjustment",
            "Quantitative SIPS vs\nPS_scz with adjustment"
        )
        for(i_test in 1:length(tests_to_plot)) {
            test <- tests_to_plot[i_test]
            main <- names_tests_to_plot[i_test]
            plot(
                x = r_g_scz_sub_range,
                y = power_mat_conditional[, test],
                xlab = "correlation (r_g)",
                ylab = "Power",
                type = "l",
                main = main,
                ylim = c(0, 1)
            )
            abline(h = 0.05, col = "red")
        }
        ## now, do adjust, no adjust
        dev.off()
    }

}


## must be run from command line
## otherwise things not sucked up!
power_analysis_aim2B <- function(fileprefix, n_power_reps = 200, nGrids = 21, model_params = get_default_model_params()) {

    r_g <- model_params$r_g
    ## Aim2B - interested in relationship between Aim2A, Aim2B
    ## what predicts iqDecline?
    ## try varying
    r1 <- seq(-1, 1, length.out = nGrids)
    r2 <- seq(-1, 1, length.out = nGrids)
    ## blank all relationships - why do I do this...
    r_g[, ] <- 0
    diag(r_g) <- 1
    r_g["scz", "iq"] <- r_g["iq", "scz"] <- -0.234
    out_mat <- array(0, c(length(r1), length(r2), 6))
    c_mat <- array(0, c(length(r1), length(r2)))
    r_mat <- array(0, c(length(r1), length(r2), 2))
    r_mat_exist <- array(FALSE, c(length(r1), length(r2)))
    for(i1 in 1:length(r1)) {
        for(i2 in 1:length(r2)) {
            print(paste0(date(), ": ", i1, " and ", i2))
            r_g["iqDecline", "iq"] <- r_g["iq", "iqDecline"] <- r1[i1]
            r_g["iqDecline", "scz"] <- r_g["scz", "iqDecline"] <- r2[i2]
            model_params$r_g <- r_g
            ## remove sub from consideration
            out <- parallel::mclapply(1:n_power_reps, mc.cores = nCores, function(iRep) {
                check <- tryCatch({
                    pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
                    0
                },
                error = function(e) 1
                )
                r <- NULL
                if (check == 0) {
                    r <- analyze_pheno_for_Aim2B(pheno) #, subpheno = subpheno)
                }
                return(list(check = check, r = r))
            })
            for(iRep in 1:n_power_reps) {
                r <- out[[iRep]]$r
                check <- out[[iRep]]$check                
                if (check == 0) {
                    dimnames(out_mat)[[3]] <- names(r)
                    out_mat[i1, i2, ] <- out_mat[i1, i2, ] + as.integer(r < 0.05)
                    c_mat[i1, i2] <- c_mat[i1, i2] + 1
                    ## r1[X, Y, 1] is iqDecline and iq
                    ## r1[X, Y, 2] is iqDecline and sz
                    r_mat[i1, i2, 1:2] <- c(r1[i1], r2[i2])
                    r_mat_exist[i1, i2] <- TRUE
                }
            }
        }
    }
    r_g <- model_params$r_g
    
    ## ## Aim2B - interested in relationship between Aim2A, Aim2B
    ## ## what predicts iqDecline?
    ## ## try varying
    ## r1 <- seq(-1, 1, length.out = nGrid)
    ## r2 <- seq(-1, 1, length.out = nGrid)
    ## ## blank all relationships
    ## r_g[, ] <- 0
    ## diag(r_g) <- 1
    ## r_g["scz", "iq"] <- r_g["iq", "scz"] <- -0.234
    ## out_mat <- array(0, c(length(r1), length(r2), 6))
    ## c_mat <- array(0, c(length(r1), length(r2)))
    ## r_mat <- array(0, c(length(r1), length(r2), 2))
    ## for(i1 in 1:length(r1)) {
    ##     for(i2 in 1:length(r2)) {
    ##         print(paste0(date(), ": ", i1, " and ", i2))
    ##         r_g["iqDecline", "iq"] <- r_g["iq", "iqDecline"] <- r1[i1]
    ##         r_g["iqDecline", "scz"] <- r_g["scz", "iqDecline"] <- r2[i2]
    ##         model_params$r_g <- r_g
    ##         ## remove sub from consideration
    ##         out <- parallel::mclapply(1:n_power_reps, mc.cores = 4, function(iRep) {
    ##             check <- tryCatch({
    ##                 pheno <- simulate_full(model_params = model_params, do_checks = FALSE)$pheno
    ##                 0
    ##             },
    ##             error = function(e) 1
    ##             )
    ##             rA <- NULL
    ##             if (check == 0) {
    ##                 rA <- analyze_pheno_for_Aim2B(pheno) ##, subpheno = "sub")
    ##             }
    ##             return(list(check = check, rA = rA)) ##, rB = rB))
    ##         })
    ##         for(iRep in 1:n_power_reps) {
    ##             rA <- out[[iRep]]$rA
    ##             ##rB <- out[[iRep]]$rB                
    ##             check <- out[[iRep]]$check                
    ##             if (check == 0) {
    ##                 dimnames(out_mat)[[3]] <- names(rA)
    ##                 out_mat[i1, i2, ] <- out_mat[i1, i2, ] + as.integer(rA < 0.05)
    ##                 c_mat[i1, i2] <- c_mat[i1, i2] + 1
    ##                 r_mat[i1, i2, 1:2] <- c(r1[i1], r2[i2])
    ##             }
    ##         }
    ##     }
    ## }
    
    save(
        r_mat,
        c_mat,
        out_mat,
        r_mat_exist,        
        file = paste0(fileprefix, ".power.stuff.RData")
    )

    name_remapper <- function(x) {
        if (x == "p_iq_vs_PRS_scz") {
            return("IQ ~ PRS_scz")
        }
        if (x == "p_iq_vs_PRS_iq") {
            return("IQ ~ PRS_iq")
        }
        if (x == "p_iqDeclineBinary_vs_PRS_scz") {
            return("(IQ decline 0/1) ~ PRS_scz")
        }
        if (x == "p_iqDeclineBinary_vs_PRS_iq") {
            return("(IQ decline 0/1) ~ PRS_iq")
        }
        if (x == "p_iqDecline_vs_PRS_scz") {
            return("IQ decline ~ PRS_scz")
        }
        if (x == "p_iqDecline_vs_PRS_iq") {
            return("IQ decline ~ PRS_iq")
        }
    }
    
    ## visualize
    plots <- lapply(1:6, function(i) {
        message(i)
        grid <- data.frame(
            X = c(r_mat[, , 1]),
            Y = c(r_mat[, , 2]),
            Z = c(out_mat[, , i] / c_mat)
        )
        r <- rasterFromXYZ(grid)
        col <- colorRampPalette(c(cbPalette[3], cbPalette[2]))(99)
        at <- seq(0, 1, length.out = 100)
        at2 <- seq(0, 1, length.out = 5)
        colorkey <- list(at = at, labels = list(labels = at2, at = at2, space = "right"))
        mapTheme <- rasterTheme(region = col)
        p <- levelplot(
            r,
            par.settings = mapTheme,
            at = at,
            colorkey = colorkey,
            margin = FALSE,
            xlab = "r_g IQ decline and IQ",
            ylab = "r_g IQ decline and SCZ",
            main = name_remapper(dimnames(out_mat)[[3]][i])
        )
        print(p)
    })

    for(suffix in c("png", "pdf")) {
        image_open(filename = fileprefix, height = 8, width = 12, suffix = suffix)
        grid.arrange(
            plots[[1]],
            plots[[5]],
            plots[[3]],            
            plots[[2]],
            plots[[6]],
            plots[[4]],
            ncol = 3, nrow = 2
        )
        dev.off()
    }

    
}



make_simple_power_matrix <- function() {

    ## make output matrix
    to_out <- array(NA, c(8, 5))
    colnames(to_out) 
    ## for scz, sub
    for(i_ps in 1:2) {
        prs_colname <- c("PS_sz", "PS_iq")[i_ps]
        prs_colnameWithR <- c("PRS_scz", "PRS_iq")[i_ps] ## lol
        if (i_ps == 1) {load(file.path(results_dir, "aim2a.power.power.stuff.RData")) } 
        if (i_ps == 2) {load(file.path(results_dir, "aim2a.iq.power.power.stuff.RData")) }
        ##
        ## scz and prs
        ##
        to_out[4 * (i_ps - 1) + 1, ] <- c(
            "scz", prs_colname, NA,
            power_mat["0", "Case_SSD_vs_Control"],
            "XXX"
        )
        ##
        ## sub and prs
        ##
        to_out[4 * (i_ps - 1) + 2, ] <- c(
            "sub", prs_colname, NA,
            paste0(
                power_mat["0", "PutativeSubthreshold_vs_Control&PutativeControl"], " [r_g = 0], ",
                power_mat["1", "PutativeSubthreshold_vs_Control&PutativeControl"], " [r_g = 1]"
            ),
            "XXX"
        )
        ##
        ## baseline FSIQ and ps
        ##
        load(file = paste0(file.path(results_dir, "aim2b.power"), ".power.stuff.RData"))
        to_out[4 * (i_ps - 1) + 3, ] <- c(
            "FSIQ", prs_colname, NA,
            round(mean(out_mat[, , paste0("p_iq_vs_", prs_colnameWithR)] / c_mat, na.rm = TRUE), 2),
            "XXX"
        )
        ## r1[X, Y, 1] is iqDecline and iq
        ## r1[X, Y, 2] is iqDecline and sz
        ## take max and 0
        w1 <- abs(r_mat[, , 2]) == max(abs(r_mat[, , 2]), na.rm = TRUE)
        w2 <- r_mat[, , 1] == 0
        w3 <- r_mat[, , 2] == 0
        ## no correlation
        w <- (w3 == TRUE) & (w3 == w2)
        w[is.na(w)] <- FALSE
        m1 <- mean(out_mat[, , paste0("p_iqDeclineBinary_vs_", prs_colnameWithR)][w] / c_mat[w])
        ## max correlation
        w <- (w1 == TRUE) & (w1 == w2)
        w[is.na(w)] <- FALSE
        m2 <- mean(out_mat[, , paste0("p_iqDeclineBinary_vs_", prs_colnameWithR)][w] / c_mat[w])
        ## 
        to_out[4 * (i_ps - 1) + 4, ] <- c(
            "viq decline", prs_colname, NA,
            paste0(
                m1, " [r_g = 0], ",
                m2, " [r_g = 1]"
            ),
            "XXX"
        )
        ## 
    }

    write.csv(
        to_out,
        file = file.path(results_dir, "power.simple.csv"),
        row.names = FALSE,
        quote = TRUE
    )

    
    
    
}

