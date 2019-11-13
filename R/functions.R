cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

run_PRSice_all_phenotypes <- function(target_file, rebuild = TRUE, nCores = 1) {

    for(phenotype in names(pheno_info)) {

        base_file <- pheno_info[[phenotype]][["file"]]
        p <- pheno_info[[phenotype]]

        fam_file_in <- paste0(target_file, ".fam")
        fam_file_out <- paste0(target_file, ".fam.dummy")
        fam <- read.table(fam_file_in)
        fam[, 6] <- fam[, 5]
        write.table(fam, file = fam_file_out, row.names = FALSE, col.names =FALSE, sep = " ", quote = FALSE)

        run_PRSice(
            base_file = base_file,
            target_file = paste0(target_file, ",", fam_file_out),
            out_dir = file.path(results_dir, paste0(phenotype, "_iBBC")),
            chr = p[["chr"]],
            bp = p[["bp"]],
            snp = p[["snp"]],
            effect_allele = p[["effect_allele"]],
            non_effect_allele = p[["non_effect_allele"]],
            base_phenotype_is_binary = p[["base_phenotype_is_binary"]],
            beta_col = p[["beta"]],
            or_col = p[["or"]],
            se = p[["se"]],
            pvalue = p[["pvalue"]],
            rebuild = rebuild,
            keep_ambig = p[["keep_ambig"]],
            seed = p[["seed"]]
        )
        
    }

}

run_PRSice <- function(
    base_file,
    target_file,
    binary_target = TRUE,
    effect_allele = "A1",
    non_effect_allele = "A2",
    bp = "BP",
    chr = "CHR",
    out_dir = dir(),
    info_col = "INFO",
    info_min = "0.90",
    snp = "SNP",
    beta_col = "BETA",
    or_col = "OR",
    base_phenotype_is_binary = TRUE,
    se = "SE",
    pvalue = "P",
    prs_p_interval = 0.05,
    prs_p_lower = 0.00,
    prs_p_upper = 0.20,
    output_all_score = TRUE,
    rebuild = TRUE,
    keep_ambig = FALSE,
    score = "std",
    seed = 3
) {
    ## defaults of PRSice
    ##    prs_p_interval = 0.000050,
    ##prs_p_lower = 0.0001,
    ##prs_p_upper = 0.50
    if (rebuild == FALSE) {
        if (file.exists(file.path(out_dir, "PRSice.all.score"))) {
            return(NULL)
        }
    }
    unlink(out_dir)
    dir.create(out_dir, showWarnings = FALSE)
    args <- c(
        "--A1", effect_allele,
        "--A2", non_effect_allele,
        "--base", shQuote(base_file),
        "--target", shQuote(target_file),
        "--bp", bp,
        "--chr", chr,
        "--info-base", paste0(info_col, ",", info_min),
        "--interval", prs_p_interval,
        "--lower", prs_p_lower,
        "--upper", prs_p_upper,
        "--se", se,
        "--snp", snp,
        "--pvalue", pvalue,
        "--prsice", shQuote(PRSice_exec),
        "--out", shQuote(file.path(out_dir, "PRSice")),
        "--seed", seed,
        "--score",  score
    )
    if (keep_ambig) {
        args <- c(args, "--keep-ambig")
    }
    if (base_phenotype_is_binary) {
        args <- c(args, "--stat", or_col)        
    } else {
        args <- c(args, "--beta", "--stat", beta_col)
    }
    if (output_all_score) {
        args <- c(args, "--all-score")
    }
    ## /home/rwdavies/bin/PRSice/downloaded/PRSice_linux \
    ##     --A1 A1 \
    ##     --A2 A2 \
    ##     --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    ##     --base /hpf/largeprojects/tcagstor/scratch/rwdavies/iBBC_22Q11/external/iPSYCH-PGC_ASD_Nov2017 \
    ##     --binary-target T \
    ##     --bp BP \
    ##     --chr CHR \
    ##     --clump-kb 250 \
    ##     --clump-p 1.000000 \
    ##     --clump-r2 0.100000 \
    ##     --info-base INFO,0.9 \
    ##     --interval 0.000050 \
    ##     --lower 0.000100 \
    ##     --model add \
    ##     --out PRSice \
    ##     --pvalue P \
    ##     --se SE \
    ##     --seed 335560616 \
    ##     --snp SNP \
    ##     --stat OR \
    ##     --target /hpf/largeprojects/tcagstor/scratch/rwdavies/iBBC_22Q11/external/22q_IBBC_forRobbie_Imputed_QCed \
    ##     --thread 1 \
    ##     --upper 0.500000
    ## argh - bugfix for now, make "relative"?
    ## dammit R and dammit PRSice
    ## a <- system2(shQuote(PRSice_R), args)     ## DAMMIT R
    x <- system(
        paste0(
            shQuote(PRSice_R),
            " ",
            paste0(args, collapse = " ")
        )
    )
    if (x > 0) {
        stop("Failed run, see above log")
    }
}



gn <- function(x, text, exl = NULL) {
    i <- grep(x, text)
    if (is.null(exl) == FALSE) {
        j <- grep(exl, text[i])
        if (length(j) > 0) {
            i <- i[-j]
        }
    }
    if (length(i) == 0) {
        stop(paste0("cannot find:", x))
    }
    if (length(i) > 1) {
        stop(paste0("multiple matches:", paste0(text[i], collapse = "---")))
    }
    return(i)
}

st <- function(i, text) {
    r <- as.numeric(strsplit(text[i], " ")[[1]][1])
    return(r)
}

get_numbers_from_log <- function(log_file) {
    text <- readLines(log_file)
    to_out <- NULL
    ## base file
    i <- gn("variant\\(s\\) included", text, exl = "base file")
    to_out <- c(
        to_out,
        n_original_in_target = st(i, text),
        n_ambiguous_target_excluded = st(i - 1, text)
    )
    ## argh this
    i <- gn("observed in base file, with:", text)
    to_out <- c(
        to_out,
        n_total_in_base = st(i, text),
        n_ambiguous_base_excluded = st(i + 1, text),
        n_base_not_in_target = st(i + 2, text),
        n_mismatched_excluded = st(i + 3, text),
        n_incl_from_base = st(i + 4, text)
    )
    i <- gn("Number of variant\\(s\\) after clumping ", text)
    to_out <- c(
        to_out,
        n_after_clumping = st(i, text)
    )
    return(to_out)
}


add_all_prs_to_pheno <- function(default_pval_cutoff = "X0.100000") {
    for(phenotype in names(pheno_info)) {
        x <- pheno_info[[phenotype]]
        if ("pval_cutoff" %in% names(x)) {
            pval_cutoff <- x[["pval_cutoff"]]
        } else {
            pval_cutoff <- default_pval_cutoff
        }
        all_score_file <- file.path(
            results_dir,
            paste0(phenotype, "_iBBC"),
            "PRSice.all.score"
        )
        all_score <- read.table(all_score_file, header = TRUE)
        a <- all_score[, c("IID", pval_cutoff)]
        colnames(a)[2] <- paste0("PRS_", phenotype)
        pheno <- merge(pheno, a, by = "IID")
    }
    return(pheno)
}



match_iid_to_affy_ids <- function(iid, pheno) {
    f <- function(x) {
        if (length(x) == 0) {
            x <- NULL
        }
        return(x)
    }
    x <- sapply(
        iid,
        function(x) {
            y1 <- f(grep(x, pheno[, "affy_id_1"]))
            y2 <- f(grep(x, pheno[, "affy_id_2"]))
            y3 <- f(grep(x, pheno[, "affy_id_3"]))
            return(c(y1, y2, y3))
        }
    )
    if (sum(is.na(x)) > 0) {
        stop("bad assumption sierwoierhwer")
    }
    return(x)
}





print_log_counts <- function() {
    ## maybe just transport PRS over then?
    sapply(names(pheno_info), function(phenotype) {
        ## get number of variants, when they were removed
        ## write to table output? 
        ## plot distribution, p-value
        print(phenotype)
        log_file <- file.path(results_dir, paste0(phenotype, "_iBBC"), "PRSice.log")
        return(get_numbers_from_log(log_file))
    })
}


check_af_between_iBBC_and_neale <- function() {
    png(paste0(intersect_out_file, ".IBBC.neale.png"), height = 5, width = 5, units = "in", res = 300)
    ## check allele frequency
    af1 <- iBBC_snps2[w1, "af"]
    af2 <- out[[2]][w3, "AC"] / out[[2]][w3, "nCompleteSamples"] / 2
    ## flip?
    af2[a1 != a3] <- 1 - af2[a1 != a3]
    m <- round(seq(1, length(af1), length.out = 1e4))
    plot(af1[m], af2[m], xlim =c(0, 1), ylim = c(0, 1))
    ## perfect fit - no need to remove un-ambiguous variants for neale and iBBC
    dev.off()
}

image_open <- function(filename, suffix = "png", height = 8, width = 12, res = 300, units = "in", remove_suffix = FALSE) {
    if (remove_suffix) {
        filename <- substr(filename, 1, nchar(filename) - 4)
    }
    filename <- paste0(filename, ".", suffix)
    if (suffix == "png") {
        png(filename = filename, width = width, height = height, res = res, units = units)
    } else {
        pdf(file = filename, width = width, height = height)
    }
}
