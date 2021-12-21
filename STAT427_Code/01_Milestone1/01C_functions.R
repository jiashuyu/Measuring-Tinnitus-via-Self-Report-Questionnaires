library(matrixcalc)
library(psych)
library(lavaan)

# Purpose: Helper functions for "01C_efacfa_analysis.Rmd"

get_efa_gof <- function(efa, gof_type){
  # Purpose: Get goodness of fit of EFA model
  return(gof_lookup[["efa"]][[gof_type]](efa))
}

######################
###### EFA Goodness of Fit
######################


get_efa_sse <- function(efa){
  # Purpose: Calculate Squared difference between sample and estimated variance
  res = efa$residual
  return(2*sum((res[lower.tri(res)])^2))
}

get_efa_loglik <- function(efa){
  # Purpose: Calculate (-2/(np))*logLik assuming normality
  n.items = dim(efa$r)[2]
  ll = loglike_mvnorm(
    M=rep(0,n.items), S=efa$r
    , mu=rep(0,n.items), Sigma=efa$stat424_varEst
    , n=nobs(efa)
  )
  return(ll * (-2/(nobs(efa)*n.items)))
}

get_efa_mledist <- function(efa){
  # Purpose: Calculate -ln(det(EstVar^-1 * SampleVar)) + tr(EstVar^-1 * SampleVar)-n.items.
  # Note: Equivalent objective to normal MLE
  # fa$objective = log(trace ((FF'+U2)^{-1} R) - log(|(FF'+U2)^-1 R|) - n.items
  # n.items = dim(efa$r)[2]
  return(efa$objective)
}

get_efa_resgen <- function(efa, V=NULL){
  # Purpose: tr[((SampleVar-EstVar)V)^2]
  # Note: if V=I or V=1 then proportional to get_efa_sse
  
  if (is.null(V)){
    n.items = dim(efa$r)[2]
    V = diag(n.items)
  }
  
  sampleVar = efa$r
  estVar = efa$stat424_varEst
  return(tr(matrix.power((sampleVar-estVar)%*%V,2)))
}

get_efa_cumss <- function(efa){
  # Purpose: Get "Cumulative Variance Explained (%)
  
  # 1-sum(efa$uniquenesses)/tr(efa$r)
  return(sum(efa$Vaccounted["Proportion Var",])) 
}

get_efa_cumss2 <- function(efa){
  # Purpose: Get Cumulative Variance UnExplained (%)
  
  return(1-get_efa_cumss(efa)) 
}

######################
###### CFA Goodness of Fit
######################


get_cfa_gof <- function(cfa, gof_type="mle_lik"){
  # Purpose: Get goodness of fit of CFA model
  return(gof_lookup[["cfa"]][[gof_type]](cfa))
}

get_cfa_sse <- function(cfa){
  # Purpose: Calculate Squared difference between sample and estimated variance
  res = resid(cfa)$cov
  return(mean((res[lower.tri(res)])^2))
}

get_cfa_loglik <- function(cfa){
  # Purpose: Calculate (-2/(np))*logLik assuming normality
  # loglike_mvnorm(M=rep(0,dim(fitted(cfa)$cov)[1]), S=attr(attr(cfa, "SampleStats"), "cov")[[1]], mu=rep(0,dim(fitted(cfa)$cov)[1]), Sigma=fitted(cfa)$cov, n=nobs(cfa))
  n.items = attr(attr(cfa,"Model"), "nvar")
  return(as.numeric(logLik(cfa) * (-2/(nobs(cfa)*n.items))))
}

get_cfa_mledist <- function(cfa){
  # Purpose: Calculate -ln(det(EstVar^-1 * SampleVar)) + tr(EstVar^-1 * SampleVar)-n.items.
  # Note: Equivalent objective to normal MLE
  n.items = attr(attr(cfa,"Model"), "nvar")
  estVar = fitted(cfa)$cov
  sampleVar = attr(attr(cfa, "SampleStats"), "cov")[[1]]
  mat = solve(estVar) %*% sampleVar
  return(-log(det(mat)) + tr(mat) - n.items)
}

get_cfa_resgen <- function(cfa, V=NULL){
  # Purpose: tr[((SampleVar-EstVar)V)^2]
  # Note: if V=I or V=1 then proportional to get_efa_sse
  
  if (is.null(V)){
    n.items = dim(fitted(cfa)$cov)[1]
    V = diag(n.items)
  }
  
  sampleVar = attr(attr(cfa, "SampleStats"), "cov")[[1]]
  estVar = fitted(cfa)$cov
  return(tr(matrix.power((sampleVar-estVar)%*%V,2)))
}

get_cfa_cumss <- function(cfa){
  # Purpose: Get Cumulative Variance Explained (%)
  
  return(1-tr(lavInspect(cfa, what="est")$theta)/tr(attr(attr(cfa, "SampleStats"), "cov")[[1]]))
}


get_cfa_cumss2 <- function(cfa){
  # Purpose: Get Cumulative Variance UnExplained (%)
  
  return(1-get_cfa_cumss(cfa))
}


######################
###### 
######################


gof_lookup <- list(
  "efa"=list(
    "resid"= get_efa_sse
    , "mle_lik"=get_efa_loglik
    , "mle_dist"=get_efa_mledist
    , "resgen"=get_efa_resgen
    , "cumss2"=get_efa_cumss2
    , "cumss"=get_efa_cumss
    )
  , "cfa"=list(
    "resid"= get_cfa_sse
    , "mle_lik"=get_cfa_loglik
    , "mle_dist"=get_cfa_mledist
    , "resgen"=get_cfa_resgen
    , "cumss2"=get_cfa_cumss2
    , "cumss"=get_cfa_cumss
    )
  , "desc" = c(
    "resid"="SSE of Correlation"
    , "mle_lik"="-2LogLik/np"
    , "mle_dist"="MLE distance"
    , "resgen"="Generic Goodness of Fit"
    , "cumss2"="Total Variance UnExplained (%)"
    , "cumss"="Total Variance Explained (%)"
    )
)


######################
###### 
######################


loading_2_matrix <- function(iloading, fit_type="resid"){
  class(iloading) = "matrix"
  iloading
}

get_lims <- function(icol){
  c(min(icol), max(icol))
}