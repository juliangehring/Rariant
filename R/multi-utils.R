combineCalls <- function(l) {
    sample_names = names(l)
    if(is.null(sample_names))
        sample_names = as.character(seq_along(l))
    return(l)
}

#data = Rariant:::exampleCalls()

yesNoMaybe <- function(x, null = 0, one = 0.5) {
    is_null = ciCovers(x, null) | ciCovers(x, -null)
    is_one = ciCovers(x, one) | x$lower > one | ciCovers(x, -one) | x$upper < -one
    is_both = is_null & is_one
    is_none = !is_null & !is_one
    out = rep(NA_character_, length(x))
    lvs = c("absent", "present", "dontknow", "inbetween")
    out[is_null] = 1L
    out[is_one] = 2L
    out[is_both] = 3L
    out[is_none] = 4L
    out = factor(out)
    levels(out) = lvs
    return(out)
}

eventFillScale <- function() {
    values = c(somatic = "#fc8d62", hetero = "#e78ac3", undecided = "#8da0cb", powerless = "#a6d854", none = "white")
    return(scale_fill_manual(values = values))
}

verdictColorScale <- function() {
    values = c(present = "#fc8d62", inbetween = "#8da0cb", dontknow = "#a6d854", absent = "lightgray")
    return(scale_color_manual(values = values))
}

updateCalls <- function(x, ...) {
    x$verdict = yesNoMaybe(x, ...)
    x$loc = factor(seq_along(x))
    return(x)
}

filterCalls <- function(x, ..., minCount) {
    idx = sapply(x, function(y, ...) {y %in% subset(y, ...)}, ...)
    if(!missing(minCount))
        idx = rowSums(idx) >= minCount
    return(idx)
}
