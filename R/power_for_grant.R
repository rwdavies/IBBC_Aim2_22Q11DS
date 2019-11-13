n <- 650
h2 <- 0.05
nRep <- 1000
## fuck it, do it again
sd <- 15
## 
ps <- sapply(1:nRep, function(iRep) {
    ##
    ps <- rnorm(n = n, mean = 0, sd = sd * sqrt(h2))
    not_ps <- rnorm(n = n, mean = 0, sd = sd * sqrt(1 - h2))
    ## normal sd is 15
    iq <- 70 + ps + not_ps
    id <- (iq < 70)
    ## want to see association
    return(coefficients(summary(glm(id ~ ps, family = "binomial")))[2, 4])
})
## 80% for binary
sum(ps < 0.05) / nRep

## can we assess differences in variance?
n1 <- 650
n2 <- 650
ps <- sapply(1:nRep, function(iRep) {
    ps1 <- rnorm(n = n, mean = 0, sd = 15 * sqrt(0.05))
    ps2 <- rnorm(n = n, mean = 0, sd = 15 * sqrt(0.05 * 0.5))
    return(var.test(ps1, ps2)$p.value)
})
sum(ps < 0.05) / nRep


## Objective 2

## Objective 3


