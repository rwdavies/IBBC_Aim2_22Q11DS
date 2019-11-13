## iBBC is hg19 e.g. rs2879914 is 51,186,143

## grove_2017_ASD = c(
##     file = file.path(external_dir, "iPSYCH-PGC_ASD_Nov2017"),
##     pretty_name = "Grove 2017 ASD"
## ),
## sniekers_2017_intelligence = c(
##     file = file.path(external_dir, "sniekers_2017_sumstats.txt.gz"),
##     pretty_name = "Sniekers 2017 intelligence"
## ),
## lee_2018_EA = c(
##     file = file.path(external_dir, "GWAS_CP_all.txt.gz"),
##     pretty_name = "Lee 2018 EA"
## ),
## lee_2018_CP = c(
##     file = file.path(external_dir, "GWAS_EA_excl23andMe.txt.gz"),
##     pretty_name = "Lee 2018 CP"
## )


## no betas
##    davies_2018_gcf = c(
##        file = file.path(external_dir, "Davies_NC_2018", "Davies2018_UKB_VNR_summary_results_29052018.txt"),
##        pretty_name = "Davies 2018 General cognitive function"
##    )


get_pheno_list <- function() {

    default_params <- list(
        binary_target = TRUE,
        effect_allele = "A1",
        non_effect_allele = "A2",
        bp = "BP",
        chr = "CHR",
        out_dir = dir(),
        info = "INFO",
        info_min = "0.90",
        snp = "SNP", ## means rsid
        beta = "BETA",
        or = "OR",
        base_phenotype_is_binary = TRUE,
        se = "SE",
        pvalue = "P",
        prs_p_interval = 0.10,
        prs_p_lower = 0.10,
        prs_p_upper = 0.50,
        output_all_score = TRUE,
        rebuild = TRUE,
        keep_ambig = FALSE
    )

    PGC_2014_SCZ <- default_params
    PGC_2014_SCZ["file"] <- file.path(external_dir, "ckqny.scz2snpres.gz")
    PGC_2014_SCZ["pretty_name"] <- "PGC 2014 SCZ"
    PGC_2014_SCZ["pval_cutoff"] <- "X0.050000"
    PGC_2014_SCZ["snp"] <- "snpid"
    PGC_2014_SCZ["chr"] <- "hg19chrc"
    PGC_2014_SCZ["bp"] <- "bp"
    PGC_2014_SCZ["effect_allele"] <- "a1"
    PGC_2014_SCZ["non_effect_allele"] <- "a2"
    PGC_2014_SCZ["info"] <- "info"
    PGC_2014_SCZ["base_phenotype_is_binary"] <- TRUE    
    PGC_2014_SCZ["or"] <- "or"
    PGC_2014_SCZ["pvalue"] <- "p"
    PGC_2014_SCZ["seed"] <- 4428 ## GE volume-2 Nov 21 2018
    #PGC_2014_SCZ["keep_ambig"] <- TRUE

    PGC_2018_SCZ <- default_params
    PGC_2018_SCZ["file"] <- file.path(external_dir, "clozuk_pgc2.meta.sumstats.reformatted.txt.gz")
    PGC_2018_SCZ["pretty_name"] <- "PGC 2018 SCZ"
    ## "As in the PGC study, we found the best P-value threshold for discrimination to be 0.05 "
    PGC_2018_SCZ["pval_cutoff"] <- "X0.050000"
    ## SNP	Freq.A1	CHR	BP	A1	A2	OR	SE	P	variant	rsid
    PGC_2018_SCZ["snp"] <- "rsid"
    PGC_2018_SCZ["chr"] <- "CHR"
    PGC_2018_SCZ["bp"] <- "BP"
    PGC_2018_SCZ["info"] <- NA
    PGC_2018_SCZ["effect_allele"] <- "A1"
    PGC_2018_SCZ["non_effect_allele"] <- "A2"
    PGC_2018_SCZ["base_phenotype_is_binary"] <- TRUE    
    PGC_2018_SCZ["or"] <- "OR"
    PGC_2018_SCZ["pvalue"] <- "P"
    PGC_2018_SCZ["seed"] <- 0548 ## GE volume last 4 numbers, Jan 7 2019
    
    Neale_XXXX_fluid <- default_params
    Neale_XXXX_fluid["file"] <- file.path(external_dir, "fluid_intelligence.neale.tsv.gz")
    Neale_XXXX_fluid["pretty_name"] <- "UKBB Fluid intelligence (Neale et al)"
    Neale_XXXX_fluid["pval_cutoff"] <- "X0.100000"    ## from similar Davies et al 2018
    Neale_XXXX_fluid["snp"] <- "rsid"
    Neale_XXXX_fluid["chr"] <- "chr"
    Neale_XXXX_fluid["bp"] <- "bp"
    Neale_XXXX_fluid["effect_allele"] <- "effect_allele"
    Neale_XXXX_fluid["non_effect_allele"] <- "non_effect_allele"
    Neale_XXXX_fluid["base_phenotype_is_binary"] <- FALSE
    Neale_XXXX_fluid["beta"] <- "beta"
    Neale_XXXX_fluid["se"] <- "se"
    Neale_XXXX_fluid["pvalue"] <- "pval"
    Neale_XXXX_fluid["info"] <- "info"
    Neale_XXXX_fluid["seed"] <- 3095 ## GE volume-2 Nov 2018
    ##Neale_XXXX_fluid["keep_ambig"] <- TRUE

    Davies_2018_G <- default_params

    pheno_info <- list(
        PGC_2014_SCZ = PGC_2014_SCZ,
        Neale_XXXX_fluid = Neale_XXXX_fluid,
        PGC_2018_SCZ = PGC_2018_SCZ        
    )

    
    return(pheno_info)
}

pheno_info <- get_pheno_list()
