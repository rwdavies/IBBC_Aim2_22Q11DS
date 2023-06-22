## make a plot like what was done before
## but for arbitrary percent variance explained, thresholds, etc


## checks
prs_explains_var <- 0.05
a <- sqrt(prs_explains_var)
b <- sqrt(1 - a ** 2)
var(a * rnorm(10000) + b * rnorm(10000))


get_ppv <- function(
    pop_prev = 0.01,
    prs_explains_var = 0.05,
    ile_thresh = 0.10
) {
    ## suppose there is a pop_prev population prevalence of k e.g. 1% this gives
    k <- qnorm(p = 1 - pop_prev)
    ## note that this should hold if our "pop" is a sub-population, of vanishingly small size,
    ##   that has a bump to their liability, e.g. they all have values of 1 higher
    ##   this effectively brings the threshold down
    ## so suppose Y = a * Y_1 + b * Y_2 is the total risk
    ## where Y_1 is the PRS which is N(0, 1)
    ## and Y_2 is the remaining risk also N(0, 1)
    ## now we need a and b such that
    ## a is the sqrt of prs_explains_var
    ## and b is sqrt of 1 - prs_explains_var
    a <- sqrt(prs_explains_var)
    b <- sqrt(1 - a ** 2)
    ## now suppose we know someone has a PRSile of some value, say ile_thresh e.g. 0.10
    ## so this gives Y_1 = qnorm(ile_thresh) = c
    c <- qnorm(ile_thresh)
    ## then we want to know the ppv i.e. what percent of them will get the thing given the condition
    ## Y > k | Y_1 = c
    ## a Y_1 + b Y_2 > k | Y_1 = c
    ## a * c + b Y_2 > k
    ## Y_2 > (k - a * c) / b
    ## which is
    ppv <- 1 - pnorm(q = (k - a * c) / b)
    ## note we expect the PS to be neutral when
    neutral_ile <- pnorm(k * (1 - b) / a)
    return(c(ppv, neutral_ile))
}

## choose some cutoffs

make_prevalence_focused_plot <- function(
    pop_prev,
    col = "black",
    x = NA
) {
    ## specify which of these we are interested in
    prs_explains_vars <- c(0.025, 0.05, 0.075, 0.10, 0.15, 0.20)
    ## similarly for the psiles
    if (is.na(x[1])) {
        xx <- log(50, 2) + -4:-1
        x <- 2 ** (xx)
        x <- c(x, 50, 100 - x[length(x):1]) / 100
        at_x <- log(50, 2) + -4:4
        is_log2 <- TRUE
    } else {
        at_x <- x ## hope this is OK!
        is_log2 <- FALSE
    }

    ##
    xlim <- range(at_x)
    ylim <- c(0, 1.0)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, col = "white", xlab = "PSile value", ylab = "PPV", main = paste0("Prevalence = ", pop_prev))
    axis(1, labels = FALSE, tick = TRUE, at = at_x)
    ## axis(1, at = plot_x, labels = round(iles * 100, 2), srt = 45)
    text(
        x = at_x,
        y = ylim[1] - 0.06,
        labels = round(x * 100, 2),
        xpd = NA,
        srt = 35
    )
    axis(2)
    lwds <- seq(1, 3, length.out = length(prs_explains_vars))
    ##
    for(i_prs_explains_vars in 1:length(prs_explains_vars)) {
        prs_explains_var <- prs_explains_vars[i_prs_explains_vars]
        vals <- sapply(x, function(ile_thresh) {
            out <- get_ppv(
                pop_prev = pop_prev,
                prs_explains_var = prs_explains_var,
                ile_thresh = ile_thresh
            )
        })
        neutral <- vals[2, 1]
        vals <- vals[1, ]
        lines(x = at_x, vals, lwd = lwds[i_prs_explains_vars], col = col)
        if (is_log2) {
            if (neutral < 0.5) {
                p <- log(100 * neutral, 2)
            } else {
                p <- 2 * log(50, 2) - log(100 * (1 - neutral), 2)
            }
        } else {
            ## plot directly
            p <- neutral
        }
        points(x = p, y = pop_prev, pch = 4, col = col)
        ## check neutral
    }
    legend("topleft", lwd = lwds, legend = prs_explains_vars)
    abline(h = pop_prev, lwd = 2, col = col)
}


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf("~/Downloads/for_jacob.pdf", height = 12, width = 8)
par(mfrow = c(3, 2))
pop_prevs <- c(0.02, 0.05, 0.10, 0.25, 0.40, 0.60)
for(i_pop_prev in 1:length(pop_prevs)) {
    make_prevalence_focused_plot(
        pop_prev = pop_prevs[i_pop_prev],
        col = cbPalette[i_pop_prev]
    )
}
dev.off()


