## compare UKBB, Davies results
library("data.table")
setwd("~/IBBC/external/")

## 
davies <- fread("Davies_NC_2018_OPEN_DATASET/Davies2018_OPEN_DATASET_summary_results.txt", header = TRUE, data.table = FALSE)
davies[, "Effect_allele"] <- toupper(davies[, "Effect_allele"])
davies[, "Other_allele"] <- toupper(davies[, "Other_allele"])
## beta, se will require allele frequency

## 
fluid <- fread(cmd = "gunzip -c fluid_intelligence.neale.tsv.gz", header = TRUE, data.table = FALSE, sep = "\t")

## intersect
both <- (intersect(davies[, "Markername"], fluid[, "rsid"]))
a1 <- match(both, davies[, "Markername"])
a2 <- match(both, fluid[, "rsid"])
## now, at matching, make agree
swap <- davies[a1, "Effect_allele"] != fluid[a2, "effect_allele"]
davies[a1[swap], c("Effect_allele", "Other_allele")] <- davies[a1[swap], c("Other_allele", "Effect_allele")]
davies[a1[swap], "Zscore"] <- -davies[a1[swap], "Zscore"]
## re-name
## colnames(davies)[colnames(davies) == "Effect_allele"] <- "allele1"
## colnames(davies)[colnames(davies) == "Other_allele"] <- "allele2"
davies$variant <- paste0(
    davies[, "Chromosome"], ":",
    davies[, "Position"], ":",
    toupper(davies[, "Other_allele"]), ":",  ## follow previous
    toupper(davies[, "Effect_allele"])
)


## un-pack assuming 

## check position on rsid

## what is intersection?
## compare p-values?
both <- (intersect(davies$variant, fluid$variant))
message(paste0("The Davies dataset has ", nrow(davies), " SNPs"))
message(paste0("The Neale fluid intelligence dataset has ", nrow(fluid), " SNPs"))
message(paste0("They intersect at ", length(both), " SNPs"))


## directionalize them?
b1 <- paste0(toupper(davies[a1, "Effect_allele"]), "_", toupper(davies[a1, "Other_allele"]))
b2 <- paste0(toupper(fluid[a2, "effect_allele"]), "_", toupper(fluid[a2, "non_effect_allele"]))
table(b1, b2)

## check direction of effect
p1 <- -log10(davies[a1, "P-value"])
p2 <- -log10(fluid[a2, "pval"])
z1 <- davies[a1, "Zscore"]
z2 <- fluid[a2, "beta"] / fluid[a2, "se"]
w <- (p1 > 3) | (p2 > 3)
png("~/IBBC/compare_davies_neale.png", height = 10, width = 20, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(p1[w], p2[w], xlab = "Davies", ylab = "UKBB / Neale")
plot(z1[w], z2[w], xlab = "Davies", ylab = "UKBB / Neale", xlim = c(-12, 12), ylim = c(-12, 12))
dev.off()

a <- davies
colnames(a)[3] <- "rsid"
both <- merge(
    a,
    fluid,
    by = "rsid"
)
## check p-values?


fluid[fluid[, "rsid"] == "rs2352974", ]
davies[davies[, "Markername"] == "rs2352974", ]

davies[davies[, "Markername"] == "rs2352974", ]
        Chromosome Position Markername Effect_allele Other_allele   P-value
7970047          3 49890613  rs2352974             T            C 6.666e-20
        Zscore        variant
7970047 -9.133 3:49890613:T:C
> fluid[fluid[, "rsid"] == "rs2352974", ]
               variant      rsid nCompleteSamples     AC    ytx       beta
8895075 3:49890613:C:T rs2352974           108818 106725 655779 -0.0745919
                se    tstat        pval chr      pos non_effect_allele
8895075 0.00903203 -8.25859 1.49029e-16   3 49890613                 C
        effect_allele     info        af
8895075             T 0.970399 0.0131451
> 


quit()

## new data
Meta-analysis of schizophrenia GWAS data from samples of European ancestry (N=105,318; 40,675 cases and 64,643 controls)

scz <- fread(cmd = paste0("gunzip -c ", file.path(external_dir, "ckqny.scz2snpres.gz")), data.table = FALSE)
scz[order(scz[, "p"])[1:5], ]

scz[order(scz[, "p"])[1:5], ]
        hg19chrc       snpid a1 a2       bp  info      or     se         p ngt
3475715     chr6 rs115329265  A  G 28712247 0.833 1.21288 0.0164 3.861e-32   0
3475555     chr6 rs114507210  T  G 28684183 0.849 0.83036 0.0161 6.623e-31   0
3475599     chr6 rs116594362  A  C 28689672 0.849 1.20310 0.0161 1.166e-30   0
3477364     chr6 rs145283874  T  C 29144532 0.730 0.79461 0.0200 1.774e-30   0
3475658     chr6 rs115923370  T  G 28700352 0.845 0.83152 0.0161 2.274e-30   0
