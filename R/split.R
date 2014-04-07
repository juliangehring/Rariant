splitSampleBinom <- function(x, n) {

    n1 = round(n/2 + 0.15*n^(3/4))
    n2 = n - n1

    stopifnot(all(n1 > 0) & all(n2 > 0))

    p = x/n

    x1 = rbinom(length(x), n1, p)
    x2 = x - x1

    ps = 0.5 * (x1/n1 + x2/n2)    

    return(ps)
}
