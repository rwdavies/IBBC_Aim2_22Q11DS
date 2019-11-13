library("data.table")
a <- fread(cmd = paste0("gunzip -c '/Users/robert davies/IBBC/external/clozuk_pgc2.meta.sumstats.txt.gz'"))
b <- fread(cmd = paste0("gunzip -c '/Users/robert davies/IBBC/external/ckqny.scz2snpres.gz'"))

## compare top SNPs from b in a

a2 <- a[a[, "CHR"] == 6, ]
y <- sapply(b[order(b[, "p"])[1:10], "snpid"], function(x) grep(x, a2[, 1]))

## times 10 to the minus 44
a[order(a[, "P"])[1:10], ]
## versus times 10 to the minus 32
b[order(b[, "p"])[1:10], ]

