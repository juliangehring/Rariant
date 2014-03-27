ciAdjustLevel <- function(eta0, conf_level) {

    cl = 1 - (1 - eta0) * (1 - conf_level)

    return(cl)
}
