scoreTest <- function(x1, n1, x2, n2) {
 
    p1hat = x1 / n1
    p2hat = x2 / n2
    dhat = p1hat - p2hat
    
    p0hat = (x1 + x2) / (n1 + n2)
    t = dhat / sqrt( p0hat * (1 - p0hat) * (1/n1 + 1/n2) )

    pval = 2 * pnorm(abs(-t))

    res = data.frame(dhat = dhat, p1 = p1hat, p2 = p2hat, tval = t, pval = pval)

    return(res)
}


nmTestOld <- function(x1, n1, x2, n2, delta = 0) {

    stopifnot(length(unique(sapply(list(x1, x2, n1, n2), length))) == 1)

    p1x = x1/n1
    p1y = x2/n2
    mu = p1x - p1y - delta
    idx = (mu != 0) ## TODO: subset
    fmdiff = numeric(length(p1x))
    p1x = p1x[idx]
    n1 = n1[idx]
    p1y = p1y[idx]
    n2 = n2[idx]
    mu = mu[idx]
    ## calc
    t = n2/n1
    a = 1 + t
    b = -(1 + t + p1x + t * p1y + delta * (t + 2))
    c = delta * delta + delta * (2 * p1x + t + 1) + p1x + t * p1y
    d = -p1x * delta * (1 + delta)
    v = (b/a/3)^3 - b * c/(6 * a * a) + d/a/2
    s = sqrt((b/a/3)^2 - c/a/3)
    u = ifelse(v > 0, s, -s)
    w = (pi + acos(v/u^3))/3
    p1d = 2 * u * cos(w) - b/a/3
    p2d = p1d - delta
    n0 = n1 + n2
    var = (p1d * (1 - p1d) / n1 + p2d * (1 - p2d) / n2)
    norm = n0 / (n0 - 1)
    fmdiff[idx] = mu / sqrt(var * norm)

    ## pvals
    pval = 2 * pnorm(abs(-fmdiff))
    res = data.frame(dhat = mu, p1 = p1d, p2 = p2d, tval = fmdiff, pval = pval)
    
    return(res)
}


feTest <- function(x1, n1, x2, n2, ...) {
    
    relErr = 1 + 1e-07
    pval = mapply(function(x1, x2, n1, n2) {
        k = x1 + x2
        lo = max(0, k - n2)
        hi = min(k, n1)
        support = lo:hi
        d = dhyper(support, n1, n2, k, log = TRUE)
        d = exp(d - max(d))
        d = d/sum(d)
        sum(d[d <= d[x1 - lo + 1] * relErr])
    }, x1, x2, n1, n2)   
    ##pval = VariantTools:::fisher_p(x1, x2, n1, n2, ...)
    
    return(pval)
}


nmTest <- function(x1, n1, x2, n2, delta = 0) {
  
    c1 = x1
    c2 = x2
    s1 = n1
    s2 = n2
    r1 = c1/s1
    r2 = c2/s2
    c = c1 + c2
    s = s1 + s2
    rd = r1 - r2

    l3 = s
    l2 = (s1 + 2 * s2) * rd - s - c
    l1 = (s2 * rd  - s - 2 * c2) * rd + c
    l0 = c2 * rd * (1 - rd)
    
    q = l2^3 / (3 * l3)^3 - l1 * l2 / (6 * l3^2) + l0 / (2 * l3)
    p = sqrt( l2^2 / (3 * l3)^2 - l1 / (3 * l3) ) # ?
    p = ifelse(q > 0, p, -p)
    a = 1/3 * (pi * acos( q / p^3 ) ) #

    r2d = 2 * p * cos(a) - l2 / (3 * l3) #
    r1d = r2d + rd #

    mu = r1 - r2
    v = sqrt( (r1d * (1-r1d) / s1) + r2d * (1-r2d) / s2 )
    nf = sqrt(1 - 1 / (s1 + s2))
    
    res = data.frame(p1hat = r1d, p2hat = r2d, tval = mu/v, tadj = mu/v*nf)
    
    return(res)
}
