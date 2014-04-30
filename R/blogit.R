blogit <- function(x, limit = 25, ...) {
    limit = abs(limit)
    x = logit(x, ...)
    x[x < -limit] = -limit
    x[x > limit] = limit
    return(x)
}
