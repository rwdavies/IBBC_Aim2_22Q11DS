## make a plot like what was done before
## but for arbitrary percent variance explained, thresholds, etc


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
    b <- sqrt(1 - prs_explains_var)
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
    ppv
}

## choose some cutoffs

make_prevalence_focused_plot <- function(
    pop_prev,
    col = "black"
) {
    ## specify which of these we are interested in
    prs_explains_vars <- c(0.025, 0.05, 0.075, 0.10, 0.15, 0.20)
    ## similarly for the psiles
    xx <- log(50, 2) + -4:-1
    x <- 2 ** (xx)
    iles <- c(x, 50, 100 - x[length(x):1]) / 100
    plot_x <- log(50, 2) + -4:4
    ##
    xlim <- range(plot_x)
    ylim <- c(0, 0.5)
    plot(x = 0, y = 0, xlim = xlim, ylim = ylim, axes = FALSE, col = "white", xlab = "PSile value", ylab = "PPV", main = paste0("Prevalence = ", pop_prev))
    axis(1, labels = FALSE, tick = TRUE, at = plot_x)
    ## axis(1, at = plot_x, labels = round(iles * 100, 2), srt = 45)
    text(
        x = plot_x,
        y = ylim[1] - 0.06,
        labels = round(iles * 100, 2),
        xpd = NA,
        srt = 35
    )
    axis(2)
    lwds <- seq(1, 3, length.out = length(prs_explains_vars))
    ##
    for(i_prs_explains_vars in 1:length(prs_explains_vars)) {
        prs_explains_var <- prs_explains_vars[i_prs_explains_vars]
        vals <- sapply(iles, function(ile_thresh) {
            get_ppv(
                pop_prev = pop_prev,
                prs_explains_var = prs_explains_var,
                ile_thresh = ile_thresh
            )
        })
        print(vals)
        lines(x = plot_x, vals, lwd = lwds[i_prs_explains_vars], col = col)
    }
    legend("topleft", lwd = lwds, legend = prs_explains_vars)
}


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
pdf("~/Downloads/for_jacob.pdf", height = 10, width = 5)
par(mfrow = c(2, 1))
make_prevalence_focused_plot(
    pop_prev = 0.02,
    col = cbPalette[2]
)
make_prevalence_focused_plot(
    pop_prev = 0.2,
    col = cbPalette[4]
)
dev.off()
