#' @export
add_outlier <- function(x,
                        outlier.type,
                        outlier.date,
                        outlier.name = sprintf("%s (%s)", outlier.type, outlier.date),
                        outliers.coef = 0){
  UseMethod("add_outlier", x)
}
#' @export
add_outlier.JD3_REGARIMA_SPEC <- function(x,
                                          outlier.type,
                                          outlier.date,
                                          outlier.name = sprintf("%s (%s)", outlier.type, outlier.date),
                                          outliers.coef = 0){
  # data.frame to recycle arguments
  new_out <- data.frame(outlier.type, outlier.date, outlier.name, outliers.coef)
  new_out <- as.list(new_out)
  new_out <- mapply(rjd3modelling::createOutlier,
                    as.list(new_out)[[1]],
                    as.list(new_out)[[2]],
                    as.list(new_out)[[3]],
                    as.list(new_out)[[4]],
                    SIMPLIFY = FALSE)
  names(new_out) <- NULL
  x$regression$outliers <- c(x$regression$outliers,
                             new_out)
  all_out = t(simplify2array(x$regression$outliers)[c("pos","code"),])
  dupl_out <- duplicated(all_out,fromLast = TRUE)
  if(any(dupl_out)){
    warning("Duplicated outliers removed: last outliers kept")
    x$regression$outliers <- x$regression$outliers[!dupl_out]
  }
  x
}
