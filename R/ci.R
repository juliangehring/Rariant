acCi <- function(x1, n1, x2, n2, conf_level = 0.95, clip = TRUE, split = FALSE) {

    ## estimate with the unadjusted counts
    p1hat = x1 / n1
    p2hat = x2 / n2
    dhat = p1hat - p2hat

    z = qnorm(1 - (1-conf_level)/2) ## two-sided
    
    ## adjust the counts only for the CIs
    zeta = z^2 / 4 ## t/4
    n1 = n1 + zeta*2
    n2 = n2 + zeta*2
    x1 = x1 + zeta
    x2 = x2 + zeta

    if(split) {
        p1t = splitSampleBinom(x1, n1)
        p2t = splitSampleBinom(x2, n2)
    } else {
        p1t = x1 / n1
        p2t = x2 / n2
    }

    ## shrinkage point estimate
    dt = p1t - p2t

    zi = z * sqrt( p1t * (1 - p1t) / n1 + p2t * (1 - p2t) /n2 )

    cil = dt - zi
    ciu = dt + zi

    if(clip) {
        cil = pmax(cil, -1)
        ciu = pmin(ciu, 1)
    }

    ## 'is.na' causes special cases here
    idx0 = !(is.finite(dhat))
    if(any(idx0)) {
        cil[idx0] = ciu[idx0] = NaN
    }

    res = data.frame(p1 = p1hat, p2 = p2hat, d = dhat, ds = dt, lower = cil, upper = ciu)
  
    return(res)
}
  

nhsCi <- function(x1, n1, x2, n2, conf_level = 0.95) {

    quant = qnorm(1 - (1-conf_level)/2) ## two-sided
    X = x1
    Y = x2
    m <- n1
    n <- n2
    pX <- x1/m
    pY <- x2/n
    estimate <- pX - pY
    CIX <- .wilsonScore(n = m, Y = X, quant = quant)
    lX <- CIX[ ,"l"]
    uX <- CIX[ ,"u"]
    CIY <- .wilsonScore(n = n, Y = Y, quant = quant)
    lY <- CIY[ ,"l"]
    uY <- CIY[ ,"u"]

    cil = estimate - quant * sqrt((lX * (1 - lX)/m) + (uY * (1 - uY)/n))
    ciu = estimate + quant * sqrt((uX * (1 - uX)/m) + (lY * (1 - lY)/n))

    res = data.frame(p1 = pX, p2 = pY, d = estimate, lower = cil, upper = ciu)
    
    return(res)
}

  
.wilsonScore <- function(n, Y, quant) {

    t = Y/n
    est.int = (Y + (quant^2)/2)/(n + (quant)^2)
    w.se = ((quant) * sqrt(n * t * (1 - t) + (quant^2)/4))/(n + quant^2)
    KI = cbind(l = est.int - w.se, u = est.int + w.se)

    return(KI)
}


acCi4 <- function(x1, n1, x2, n2, conf_level = 0.95, clip = TRUE) {

    ## estimate with the unadjusted counts
    p1hat = x1 / n1
    p2hat = x2 / n2
    dhat = p1hat - p2hat

    ## adjust the counts only for the CIs
    adjust = 1
    n1 = n1 + 2*adjust
    n2 = n2 + 2*adjust
    x1 = x1 + adjust
    x2 = x2 + adjust

    p1t = x1 / n1
    p2t = x2 / n2
    dt = p1t - p2t

    z = abs(qnorm((1 - conf_level)/2)) ## two-sided
    zi = z * sqrt( p1t * (1 - p1t) / n1 + p2t * (1 - p2t) /n2 )

    cil = dt - zi
    ciu = dt + zi

    if(clip) {
        cil = pmax(cil, -1)
        ciu = pmin(ciu, 1)
    }

    ## 'is.na' causes special cases here
    idx0 = !(is.finite(dhat))
    if(any(idx0)) {
        cil[idx0] = ciu[idx0] = NaN
    }

    res = data.frame(d = dhat, p1 = p1hat, p2 = p2hat, lower = cil, upper = ciu)
  
    return(res)
}
