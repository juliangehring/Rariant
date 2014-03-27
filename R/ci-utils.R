ciOutside <- function(x, delta = 0) {
    idx = x$lower > delta | x$upper < delta
    return(idx)
}


ciCovers <- function(x, delta = 0) {
    idx = x$lower <= delta & x$upper >= delta
    return(idx)
}


ciOverlap <- function(x, y) {
    idx = x$lower <= y$upper & y$lower <= x$upper
    return(idx)
}


ciWidth <- function(x) {
    width = x$upper - x$lower
    return(width)
}
