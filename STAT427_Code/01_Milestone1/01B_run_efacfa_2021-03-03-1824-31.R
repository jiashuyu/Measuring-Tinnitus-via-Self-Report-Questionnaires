rm(list=ls())
library(readxl)
# library(mirt)
library(psych)
library(lavaan)

set.seed(71)

###################
## Purpose:
## For each survey (THI,TFI,TPFQ)
## Run EFA using 1 to 8 factors
## Run CFA using existing subscales
## CFA is done with two different techniques
###################

###################
###### Get Data
###################

helper_gettab = function(tab, filename, filepath=dir_path){
  ifilepath = paste(filepath, filename, sep="/")
  read_excel(path=ifilepath, sheet=tab, na="NR")
}

dir_path = "../data" # Modify to your specific machine
tabs = c("THI", "TFI", "TPFQ")

data_all = read.csv(paste(dir_path, "merged.csv", sep="/"))
icols = c(
  which(grepl(paste0("^BAI_"), colnames(data_all)))
  , which(grepl(paste0("^BDI_"), colnames(data_all)))
)
data_all = data_all[-icols]
data_all = na.omit(data_all)


data_tabs = lapply(
  tabs
  , function(itab){
    icols = grepl(paste0("^", itab, "_"), colnames(data_all))
    data_all[icols]
  }
)
names(data_tabs) = tabs


lookups_all = lapply(
  tabs
  , function(itab) helper_gettab(itab, filename="subscore-lookup.xlsx")
)
names(lookups_all) = tabs
lookups_all[["TPFQ"]] = lookups_all[["TPFQ"]][,-2]


###################
###### Functions
###################


lookup_to_str <- function(ilookup, fixed=FALSE){
  fixed_str = ifelse(fixed,"1*", "")
  return(paste(
    sapply(1:(ncol(ilookup)), function(i){
      items = which(ilookup[,i]==1)
      istr = paste0(
        "F", i, " =~ ", fixed_str
        , paste(colnames(idata)[items]
                , collapse=paste0(" + ", fixed_str)))
      istr
    })
    , collapse="\n"))
}


###################
###### Simulate
###################


nfactor=8
nrand = 100
results = list()
for (itab in tabs){
  print(itab)
  
  idata = data_tabs[[itab]]
  ilookup = lookups_all[[itab]]
  q_legacy = ncol(ilookup)-2
  target_keys = as.matrix(ilookup[3:ncol(ilookup)]) # make.keys(nrow(ilookup), lapply(3:ncol(ilookup), function(i) which(ilookup[,i]==1)))
  spear_cor= cor(idata, use="complete.obs", method = "spearman")
  
  
  # Fit EFA
  iefa_out = list()
  for (ifactors in 1:nfactor){
    ifit = NA
    tryCatch({
      ifit  = fa(spear_cor, n.obs=nrow(idata), ifactors, fm="mle", rotate="varimax")
      ifit$stat424_varEst = ifit$r-ifit$residual + diag(ifit$uniquenesses)
      if (ifactors > 1){
        ifit$stat424_varimax = faRotate(ifit$loadings, rotate="varimax")
        ifit$stat424_oblimin = faRotate(ifit$loadings, rotate="oblimin")
      }
      if (ifactors == q_legacy){
        ifit$stat424_target = target.rot(ifit$loadings, keys=target_keys)
      }
      }, error=function(e){})
    iefa_out[[ifactors]] = ifit
  }
  
  # Fit CFA
  helper_cfa <- function(model_mat, is_reduc
                         # , spear_cor=spear_cor, nobs=nrow(idata) # Implicit
                         ){
    # Purpose: Fit CFA. Update each loop
    model_str = lookup_to_str(model_mat, is_reduc)
    cfa_fit = cfa(
      model=model_str
      , sample.cov=spear_cor # Implicit
      , sample.nobs	= nrow(idata) # Implicit nobs
    )
    attr(cfa_fit, "stat424_lookup") = model_mat
    return(cfa_fit)
  }
  all_samples = t(sapply(
    1:nrand
    , function(i) sample(nrow(ilookup),nrow(ilookup), replace=FALSE)
  ))
  cfa_list = list()
  for (itype in c("full", "reduced")){
    icfa_out = list()
    is_reduc = (itype=="reduced")
    
    # Fit Legacy CFA
    model_mat = ilookup[,2+(1:q_legacy)]
    icfa_out[["legacy"]] = helper_cfa(model_mat, is_reduc)
    
    # Legacy with just individual columns
    legacy_one = list()
    for (i in 1:q_legacy){
      # Legacy just one subscale
      model_mat = ilookup[,c(i+2)]
      legacy_one[[i]] = helper_cfa(model_mat, is_reduc)
    }
    icfa_out[["legacy_one"]] = legacy_one

    # Random Comparison
    rand_cfa = list()
    for (i in 1:nrand){
      irand = NA
      tryCatch({
        isample = all_samples[i,]
        model_mat = ilookup[isample,3:ncol(ilookup)]
        irand = helper_cfa(model_mat, is_reduc)
      } , error=function(e){})
      rand_cfa[[i]] = irand
    }
    icfa_out[["rand"]] = rand_cfa
    
    cfa_list[[itype]] = icfa_out
  }
  
  results[[itab]] = list(
    efa=iefa_out
    , cfa=cfa_list
    , spear_cor=spear_cor
    , ilookup=ilookup
  )
}

rm(list=setdiff(ls(), c(
  "tabs", "data_all", "data_tabs"
  , "lookups_all"
  , "results"
  , "lookup_to_str")))

# Save to RDATA file