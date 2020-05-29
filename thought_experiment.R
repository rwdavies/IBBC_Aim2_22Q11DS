## on sex, 1% prevalence, liability threshold model
## prs explains 5% of variance
prev <- 0.01
Z_t <- qnorm(p = 1 - prev)
var <- 1
prs_fraction_var_explained <- 0.05
s2g <- prs_fraction_var_explained * var
s2e <- (1 - prs_fraction_var_explained) * var
N <- 2e8

f <- function() {
    ## variant at 1/3000 with 25% prevalence
    has_22qDS <- runif(N) < (1 / 3000)
    g <- rnorm(n = N, mean = 0, sd = sqrt(s2g))
    Y <- g + rnorm(n = N, mean = 0, sd = sqrt(s2e))
    Y[has_22qDS] <- Y[has_22qDS] + (Z_t - qnorm(p = 0.75))
    ## Z <- Y > Z_t
    ## still basically normal with expected distribution
    mean(Y)
    var(Y)
    ## now recruit based on two strategies
    ## 1) 1000 entirely random 22q11DS
    ## 2) 75%, 25% has scz
    to_out <- lapply(1:2, function(i) {
        if (i == 1) {
            who <- sample(which(has_22qDS), 10000)
        } else if (i == 2) {
            ## first, get all of them, then sample first bit
            a <- which(has_22qDS)
            who1 <- sample(a, 7500)
            a <- setdiff(a, who1)
            who2 <- sample((Y > Z_t)[a], 2500)        
            who <- c(who1, who2)
        }
        gx <- g[who]
        Yx <- Y[who]
        Zx <- (Yx > Z_t)
        cuts <- quantile(gx, probs = seq(0, 1, length.out = 6))
        cuts[1] <- cuts[1] - 1 ## get lowest in
        out <- tapply(Zx, as.integer(cut(gx, breaks = cuts)), function(a) {
            sum(a) / length(a)
        })
        return(out)
    })
    return(to_out)
}

out <- parallel::mclapply(1:4, mc.cores = 4, function(x) f())

## plot average
outA <- out[[1]][[1]]
outB <- out[[1]][[2]]
for(i in 1:length(out)) {
    outA <- outA + out[[i]][[1]]
    outB <- outB + out[[i]][[2]]
}
outA <- outA / length(out)
outB <- outB / length(out)


    plot(out, ylim = c(0, 0.5), col = cbPalette[1 + i], type = "l")
    par(new = TRUE)
}



cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
