#' @importFrom stats pt
NULL


#' JD3 print functions
#'
#' @param x the object to print.
#' @param ... further unused parameters.
#' @name jd3_print
#' @rdname jd3_print
#' @export
print.JD3_ARIMA<-function(x, ...){
  m <- x
  if (m$var > 0 || length(m$delta)>1){
    cat(m$name, "\n\n")
    if (length(m$ar)>1) cat("AR: ", m$ar, "\n")
    if (length(m$delta)>1)cat("DIF: ", m$delta, "\n")
    if (length(m$ma)>1)cat("MA: ", m$ma, "\n")
    cat("var: ", m$var, "\n\n")
  }
  invisible(x)
}


#' @rdname jd3_print
#' @export
print.JD3_UCARIMA<-function(x,...){
  ucm <- x
  print(ucm$model)
  lapply(ucm$components, function(z){print(z)})
  invisible(x)
}

arima_node<-function(p,d,q){
  s<-paste(p,d,q,sep=',')
  return (paste0('(', s, ')'))
}


#' @rdname jd3_print
#' @export
print.JD3_SARIMA<-function(x, ...){
  m <- x
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
  invisible(x)
}


#' @rdname jd3_print
#' @export
print.JD3_SARIMA_ESTIMATION<-function(x, ...){
  tables = sarima_coef_table(x, ...)
  orders = tables$sarima_orders

  cat("SARIMA model: ",
      arima_node(orders$p, orders$d, orders$q),
      arima_node(orders$bp, orders$bd, orders$bq),
      "\n")

  cat("\nCoefficients\n")
  if(ncol(tables$coef_table) == 2){
    print(tables$coef_table)
  }else{
    print(tables$coef_table[-2])
  }
  invisible(x)
}
sarima_coef_table <- function(x, cov = NULL, ndf = NULL,...){
  m <- x

  if (! is.null(m$phi)) p<-dim(m$phi)[2]else p<-0
  if (! is.null(m$theta)) q<-dim(m$theta)[2]else q<-0
  if (! is.null(m$bphi)) bp<-dim(m$bphi)[2]else bp<-0
  if (! is.null(m$btheta)) bq<-dim(m$btheta)[2]else bq<-0
  sarima_orders = list(p = p, d = m$d, q = q, bp = bp, bd = m$bd, bq = bq)
  names<-NULL
  if (p > 0){names=c(names,paste0("phi(", 1:p, ')')) }
  if (q > 0){names=c(names,paste0("theta(", 1:q, ')')) }
  if (bp > 0){names=c(names,paste0("bphi(", 1:bp, ')')) }
  if (bq > 0){names=c(names,paste0("btheta(", 1:bq,')')) }
  if (! is.null(names)){
    all<-t(cbind(m$phi, m$theta, m$bphi, m$btheta))
    fr<-as.data.frame(all, row.names = names)
    for(i in colnames(fr)){
      fr[,i] <- unlist(fr[,i])
    }
    if(!is.null(cov) & !is.null(ndf)){
      fr$pvalue <- fr$t <- fr$stde <- NA
      stde<-sqrt(diag(cov))
      sel<-fr$type=='ESTIMATED'
      t<-fr$value[sel]/stde
      pval<-2*pt(abs(t), ndf, lower.tail = F)
      fr$stde[sel]<-stde
      fr$t[sel]<-t
      fr$pvalue[sel]<-pval
    }
  }
  list(sarima_orders = sarima_orders,
       coef_table = fr)
}


#' @rdname jd3_print
#' @export
print.JD3_SPAN<-function(x, ...){
  span <- x
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


#' @rdname jd3_print
#' @export
print.JD3_LIKELIHOOD<-function(x, ...){
  ll <- x
  cat("Number of observations: ", ll$nobs, "\n")
  cat("Number of effective observations: ", ll$neffectiveobs, "\n")
  cat("Number of parameters: ", ll$nparams, "\n\n")
  cat("Loglikelihood: ", ll$ll, "\n")
  if (ll$ll != ll$adjustedll)cat("Adjusted loglikelihood: ", ll$adjustedll, "\n\n")
  cat("Standard error of the regression (ML estimate): ", sqrt(ll$ssq/ll$neffectiveobs), "\n")
  cat("AIC: ", ll$aic, "\n")
  cat("AICC: ", ll$aicc, "\n")
  cat("BIC: ", ll$bic, "\n\n")
  invisible(x)
}


#' @rdname jd3_print
#' @export
print.JD3_REGARIMA_RSLTS<-function(x, ...){
  cat("Log-transformation:",if(x$description$log) {"yes"} else {"no"},sep=" ")
  cat("\n")
  ndf<-x$estimation$likelihood$neffectiveobs-x$estimation$likelihood$nparams+1
  print(x$description$arima, cov = x$estimation$parameters$cov,
        ndf = ndf,
        ...)
  cat("Coefficients:\n")
  cat("ARIMA:\n")
  xregs = regarima_coef_table(x, ...)
  if (!is.null(xregs)){
    cat("Regression model:\n")
    print(xregs[-2])
  }else{
    cat("No regression variables\n")
  }
  print(x$estimation$likelihood, ...)
  invisible(x)
}
regarima_coef_table <- function(x,...){
  q <- x
  if (length(q$description$variables)>0){
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
    xregs
  }else{
    NULL
  }
}
