get_and_sanitize <- function(what) {
    gsub("\\\\", "", Sys.getenv(what))
}

args <- commandArgs(trailingOnly = TRUE)

if (1 == 0) {

    setwd("~/IBBC/")
    args <- c(
        "~/IBBC/external//clozuk_pgc2.meta.sumstats.txt.gz",
        "~/IBBC/external//clozuk_pgc2.meta.sumstats.test.txt.gz"
    )
        
    ## input file
    ## output file
    
}

print(args)
input_file <- args[1]
output_file <- args[2]

library("data.table")

input <- fread(cmd = paste0("gunzip -c ", input_file), data.table = FALSE)
input$variant <- paste(input[, "CHR"], input[, "BP"], toupper(input[, "A1"]), toupper(input[, "A2"]), sep = ":")

file <- "~/IBBC/external/variants.tsv.gz"
variants <- data.table::fread(cmd = paste0("gunzip -c ", file), data.table = FALSE)

input$rsid <- paste0("rs123456789000", input[, "BP"])
m <- match(input[, "variant"], variants[, "variant"])
input[is.na(m) == FALSE, "rsid"] <- variants[m[is.na(m) == FALSE], "rsid"]

output_file_no_gz <- substr(output_file, 1, nchar(output_file) - 3)
data.table::fwrite(
    input,
    file = output_file_no_gz,
    sep = "\t",
    quote = FALSE,
    col.names = TRUE,
    row.names = FALSE
)
system(paste0("gzip -1 -f ", output_file_no_gz))

quit()

prev <- fread(cmd = paste0("gunzip -c ", "~/IBBC/external/ckqny.scz2snpres.gz"), data.table = FALSE)

input[3271, ] ## rs10000000
variants[variants[, "rsid"] == "rs10000000", ]


## match up previous and current
int <- intersect(input[, "rsid"], prev[, "snpid"])
##
m1 <- match(int, input[, "rsid"])
z1 <- log(input[m1, "OR"]) / input[m1, "SE"]
p1 <- -log10(input[m1, "P"])
##
m2 <- match(int, prev[, "snpid"])
z2 <- log(prev[m2, "or"]) / prev[m2, "se"]
p2 <- -log10(prev[m2, "p"])
##

swap <- input[m1, "A2"] != prev[m2, "a2"]
## swap alleles if not matched?
z2B <- z2
z2B[swap] <- -z2B[swap]

## plot
png("~/IBBC/compare_scz_2018_2014.png", height = 10, width = 20, units = "in", res = 300)
w <- (p1 > 3) | (p2 > 3)
par(mfrow = c(1, 2))
lim <- range(c(range(p1[w]), range(p2[w])))
plot(p1[w], p2[w], xlab = "New", ylab = "Old", xlim = lim, ylim = lim)
abline(0, 1)
b <- max(abs(range(c(range(z1[w]), range(z2[w])))))
lim <- c(-b, b)
plot(z1[w], z2B[w], xlab = "New", ylab = "Old", xlim = lim, ylim = lim)
abline(0, 1)
abline(0, -1)
dev.off()

1


input$rsid <- info[, 1]
input[input[, "rsid"] == input[, "CHR"], "rsid"] <- NA
input$bp <- info[, 2]
input$effect_allele <- info[, 3]
input$non_effect_allele <- info[, 4]



## AM HERE
## LINE UP, TEST! HOW DOES IT DO ON SCHIZOPHRENIA
int <- intersect(input[, "variant"], variants[, "variant"])
length(int)
nrow(input)
nrow(variants)



