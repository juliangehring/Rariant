yesNoMaybe <- function(x, null = 0, one = 0.5) {

    if(null >= one)
        warning("'one' should be >= than 'null'")

    is_null = ciCovers(x, null) | ciCovers(x, -null) | (x$upper < null & x$lower > -null)
    is_one = ciCovers(x, one) | x$lower > one | ciCovers(x, -one) | x$upper < -one
    is_both = is_null & is_one
    is_none = !is_null & !is_one

    out = rep(NA_character_, length(x))
    lvs = c("absent", "present", "dontknow", "inbetween")
    out[is_null] = lvs[1L]
    out[is_one] = lvs[2L]
    out[is_both] = lvs[3L]
    out[is_none] = lvs[4L]
    out = factor(out, levels = lvs)

    return(out)
}
