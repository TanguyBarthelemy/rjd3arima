
#' Title
#'
#' @param period
#' @param phi
#' @param d
#' @param theta
#' @param bphi
#' @param bd
#' @param btheta
#' @param name
#'
#' @return
#' @export
#'
#' @examples
sarima.model<-function(period, phi=NULL, d=0, theta=NULL, bphi=NULL, bd=0, btheta=NULL, name=NULL){
  return (structure(list(name=name, period=period, phi = phi, d=d, theta=theta,
                         bphi = bphi, bd = bd, btheta = btheta)), class="JD3SARIMA")
}

#' Title
#'
#' @param ar
#' @param delta
#' @param ma
#' @param var
#' @param name
#'
#' @return
#' @export
#'
#' @examples
arima.model<-function(ar=NULL, delta=NULL, ma=NULL, var=1, name=NULL){
  return (structure(list(name=name, innovationvariance=var, ar=ar, delta=delta, ma=ma), class= "JD3ARIMA"))
}

#' Title
#'
#' @param model
#' @param components
#' @param checkmodel
#'
#' @return
#' @export
#'
#' @examples
ucarima.model<-function(model, components, checkmodel=T){

# TODO: compute the model when it is missing and check the model if it is provided

if (! is(model, "JD3ARIMA") && ! is(model, "JD3SARIMA")) stop("Invalid model")
lapply(components, function(c){if (! is(c, "JD3ARIMA")) stop("Invalid component")})

return (structure
        (list(
          model=model,
          components=components),
          class= "JD3UCARIMA"))
}

