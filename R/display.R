#' @include utils.R
NULL

#' Title
#'
#' @param m
#'
#' @return
#' @export
#'
#' @examples
print.JD3ARIMA<-function(m){
  if (m$innovationvariance > 0 || length(m$delta)>1){
    cat(m$name, "\n\n")
    if (length(m$ar)>1) cat("AR: ", m$ar, "\n")
    if (length(m$delta)>1)cat("DIF: ", m$delta, "\n")
    if (length(m$ma)>1)cat("MA: ", m$ma, "\n")
    cat("var: ", m$innovationvariance, "\n\n")
  }
}


#' Title
#'
#' @param ucm
#'
#' @return
#' @export
#'
#' @examples
print.JD3UCARIMA<-function(ucm){
  print(ucm$model)
  lapply(ucm$components, function(z){print(z)})
}

arima_node<-function(p,d,q){
  s<-paste(p,d,q,sep=',')
  return (paste0('(', s, ')'))
}

#' Title
#'
#' @param m
#'
#' @return
#' @export
#'
#' @examples
print.JD3SARIMA<-function(m){
  cat("SARIMA model: ", arima_node(m$p, m$d, m$q), arima_node(m$bp, m$bd, m$bq), "\n")

  cat("\ncoefficients\n")
  if (length(m$parameters > 0)){
    names<-NULL
    e<-sqrt(diag(m$covariance))
    if (m$p > 0){names=c(names,paste("phi", 1:m$p, sep='-')) }
    if (m$bp > 0){names=c(names,paste("bphi", 1:m$bp, sep='-')) }
    if (m$q > 0){names=c(names,paste("theta", 1:m$q, sep='-')) }
    if (m$bq > 0){names=c(names,paste("btheta", 1:m$bq, sep='-')) }
    q<-data.frame(coef=m$parameters, stde=e, row.names = names)
    print(q)

    cat("\ncorrelation of the coefficients\n")
    corr<- m$covariance/e%*%t(e)
    corr<-`row.names<-`(corr, names)
    corr<-`colnames<-`(corr, names)
    print(corr)

    cat("\nscores of the coefficients\n")
    print (m$score)
  }
}


#' Title
#'
#' @param span
#'
#' @return
#' @export
#'
#' @examples
print.JD3SPAN<-function(span){
  type<-span$type
  d0<-span$d0
  d1<-span$d1
  n0<-span$n0
  n1<-span$n1

  if (type=="ALL") {x<-"All"}
  else if (type=="FROM") {x<-paste("From",d0, sep=" ")}
  else if (type=="To") {x<-paste("Until",d1, sep=" ")}
  else if (type=="BETWEEN") {x<-paste(d0,d1,sep=" - ")}
  else if (type=="FIRST") {x<-paste("All but first",n0,"periods", sep=" ")}
  else if (type=="LAST") {x<-paste("All but last",n1,"periods", sep=" ")}
  else if (type=="EXCLUDING") {x<-paste("All but first",n0,"periods and last",n1,"periods", sep=" ")}
  else {x<- "Undefined"}

  print(x)
}

#' Title
#'
#' @param ll
#'
#' @return
#' @export
#'
#' @examples
print.JD3LIKELIHOOD<-function(ll){
  cat("Number of observations: ", ll$nobs, "\n")
  cat("Number of effective observations: ", ll$neffectiveobs, "\n")
  cat("Number of parameters: ", ll$nparams, "\n\n")
  cat("Loglikelihood: ", ll$ll, "\n")
  if (ll$ll != ll$adjll)cat("Adjusted loglikelihood: ", ll$adjustedll, "\n\n")
  cat("Standard error of the regression (ML estimate): ", sqrt(ll$ssq/ll$neffectiveobs), "\n")
  cat("AIC: ", ll$aic, "\n")
  cat("AICC: ", ll$aicc, "\n")
  cat("BIC: ", ll$bic, "\n\n")
}

#' Title
#'
#' @param rslts
#'
#' @return
#' @export
#'
#' @examples
print.JD3REGARIMA_RSLTS<-function(q){

  regs<-do.call("rbind", lapply(q$description$variables, function(z){z$coeff}))
  xregs<-cbind(regs, stde=NA, t=NA, pvalue=NA)
  stde<-sqrt(diag(q$estimation$bvar))
  sel<-xregs$type=='ESTIMATED'
  t<-xregs$value[sel]/stde
  ndf<-q$estimation$likelihood$neffectiveobs-q$estimation$likelihood$nparams+1
  pval<-2*pt(abs(t), ndf, lower.tail = F)
  xregs$stde[sel]<-stde
  xregs$t[sel]<-t
  xregs$pvalue[sel]<-pval
  print(xregs[-2])
}

