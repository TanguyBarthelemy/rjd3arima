
.onLoad <- function(libname, pkgname) {
  if (! requireNamespace("rjd3modelling", quietly = T)) stop("Loading rjd3 libraries failed")
}

