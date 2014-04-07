rdSimulate <- function(pars, fun, n_sample, ...) {

    k1 = pars["k1"]
    k2 = pars["k2"]
    n1 = pars["n1"]
    n2 = pars["n2"]
    
    ks1 = rbinom(n_sample, n1, k1/n1)
    ks2 = rbinom(n_sample, n2, k2/n2)
    dhat = k1/n1 - k2/n2
    ci = fun(ks1, n1, ks2, n2, ...)
    
    idx_in = dhat >= ci$lower & dhat <= ci$upper
    cp = mean(idx_in)

    mw = mean(ci$upper - ci$lower)

    return(c(cp, mw))
}


coverageProbability <- function(pars, fun = acCi, n_sample = 1e4, min_k, ...) {
  
    if(!missing(min_k)) {
        pars$k1 = pmin(pars$k1, min_k)
        pars$k2 = pmin(pars$k2, min_k)
    }
    
    res = t(apply(pars, 1, function(x) rdSimulate(x, fun, n_sample, ...) ))
    colnames(res) = c("cp", "mw")
    
    pars$p1 = pars$k1 / pars$n1
    pars$p2 = pars$k2 / pars$n2
    pars$d = pars$p1 - pars$p2
    pars$cp = res[ ,"cp"]
    pars$aw = res[ ,"mw"]
    
    return(pars)
}
