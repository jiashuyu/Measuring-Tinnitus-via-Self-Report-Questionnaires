rm(list=ls())
library(readxl)
# library(MASS)
library(mvtnorm)

# Purpose: Fit survey A using the responses/subscales from survey B.

###################
###### Get Data
###################


helper_gettab = function(tab, filename, filepath=dir_path){
  ifilepath = paste(filepath, filename, sep="/")
  read_excel(path=ifilepath, sheet=tab, na="NR")
}

dir_path = "C:/Users/drakethrice/Box/BowersJesse/Class_STAT427_Consulting/project"
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
###### Transform Data
###################

# Data as 0-10 factors

# id_convert_10 = cbind(0L:10L,0L:10L)
# id_convert_100 = cbind(0L:100L, sapply(0L:100L, function(i) 1+floor(i/10)))
# id_convert_100[1,2] = 0
# id_convert_100[101,2] = 10
id_convert_10 = cbind(0L:10L,c(0, rep(1,5), rep(2,5)))
id_convert_100 = cbind(0L:100L,c(0, rep(1,50), rep(2,50)))

data_tabs_f <- data_tabs
for (itab in names(data_tabs_f)){
  data_tabs_f[[itab]] = as.data.frame(data_tabs_f[[itab]])
  # Convert to factor
  if (itab %in% c("TFI")){
    data_tabs_f[[itab]] = as.data.frame(apply(
      round(data_tabs_f[[itab]]), c(1,2)
      , function(x) id_convert_10[x+1,2]))
  }
  if (itab %in% c("TPFQ")){
    data_tabs_f[[itab]] = as.data.frame(apply(
      round(data_tabs_f[[itab]]), c(1,2)
      , function(x) id_convert_100[x+1,2]))
  }
  
  for (i in 1:ncol(data_tabs_f[[itab]])){
    data_tabs_f[[itab]][,i]  = factor(data_tabs_f[[itab]][,i], ordered=TRUE)
    # data_tabs_f[[itab]][,i] = as.integer(data_tabs_f[[itab]][,i])
  }
}

list2mat <- function(ilist){
  out = lapply(1:length(ilist), function(i) cbind(names(ilist)[i], ilist[[i]]))
  out = do.call(rbind, out)
  return(out)
}

# Raw scores
scores = list()
for (itab in tabs){
  iscores = (
    as.matrix(data_tabs[[itab]]) 
    %*% as.matrix(lookups_all[[itab]][,-c(1,2)]))
  colnames(iscores) = paste0(itab, "_", colnames(iscores))
  scores[[itab]] = iscores
}


###################
###### Fit Models
###################


fact2num <- function(ifactor){
  as.numeric(levels(ifactor)[ifactor])
}

# surv2surv
surv2surv = expand.grid(tabs, tabs)
colnames(surv2surv) <- c("in_tab", "out_tab")
surv2surv=surv2surv[surv2surv["in_tab"] != surv2surv["out_tab"],]
surv2surv$fit <- apply(
  surv2surv, 1
  , function(irow, df){
    tab_in = unlist(irow["in_tab"])
    tab_out = unlist(irow["out_tab"])
    my_str = paste0("cbind("
                    , paste0(colnames(df[[tab_out]]), collapse=",")
                    , ") ~ 1 + "
                    , paste0(colnames(df[[tab_in]]), collapse=" + ")
    )
    ifit = lm(
      formula(my_str)
      , data = data.frame(cbind(df[[tab_in]], df[[tab_out]]))
    )
    ifit
  }
  , df=scores
)
surv2surv$RSS <- sapply(
  surv2surv$fit
  , function(ifit) sum(resid(ifit)^2)
)
surv2surv$SS <- sapply(
  surv2surv$fit
  , function(ifit){
    y = predict(ifit) + resid(ifit)
    sum(scale(y, scale=FALSE)^2)
  }
)
surv2surv$R2 = 1-surv2surv$RSS/surv2surv$SS
surv2surv$logLik <- sapply(
  surv2surv$fit
  , function(ifit){
    sum(dmvnorm(resid(ifit), sigma=cov(resid(ifit)), log=TRUE))
  }
)
surv2surv$dof <- sapply(
  surv2surv$fit
  , function(ifit){
    dof_mean=length(ifit$coefficients)
    dof_var = dim(ifit$coefficients)[2]*(dim(ifit$coefficients)[2]+1)/2
    dof_mean+dof_var
  }
)
surv2surv$aic = -2*surv2surv$logLik+2*surv2surv$dof


# surv2score
surv2score = merge(
  tabs
  , list2mat(lapply(lookups_all, function(ilookup) colnames(ilookup)[-c(1,2)]))
  , all=TRUE
)
colnames(surv2score) = c("in_tab", "out_tab", "out_score")
surv2score=surv2score[surv2score["in_tab"] != surv2score["out_tab"],]
surv2score$fit <- apply(
  surv2score, 1
  , function(irow, df){
    tab_in = unlist(irow["in_tab"])
    tab_out = unlist(irow["out_tab"])
    score_out = unlist(irow["out_score"])
    score_out = paste0(tab_out, "_", score_out)
    my_str = paste0(score_out
                    , " ~ 1 + "
                    , paste0(colnames(df[[tab_in]]), collapse=" + ")
    )
    ifit = lm(
      formula(my_str)
      , data = data.frame(cbind(df[[tab_in]], df[[tab_out]]))
    )
    ifit
  }
  , df=scores
)
surv2score$RSS <- sapply(
  surv2score$fit
  , function(ifit) sum(resid(ifit)^2)
)
surv2score$SS <- sapply(
  surv2score$fit
  , function(ifit){
    y = predict(ifit) + resid(ifit)
    sum(scale(y, scale=FALSE)^2)
  }
)
surv2score$R2 = 1-surv2score$RSS/surv2score$SS
surv2score$dof = sapply(
  surv2score$fit
  , function(ifit)
    length(ifit$coefficients)+1
)
surv2score$AIC <- sapply(surv2score$fit, AIC)


# surv2q
surv2q = merge(
  tabs
  , list2mat(lapply(data_tabs, function(idata) colnames(idata)))
  , all=TRUE)
colnames(surv2q) = c("in_tab", "out_tab", "out_q")
surv2q=surv2q[surv2q["in_tab"] != surv2q["out_tab"],]
surv2q$fit <- apply(
  surv2q, 1
  , function(irow, df1, df2){
    tab_in = unlist(irow["in_tab"])
    tab_out = unlist(irow["out_tab"])
    q_out = unlist(irow["out_q"])
    my_str = paste0(q_out
                    , " ~ 1 + "
                    , paste0(colnames(df1[[tab_in]]), collapse=" + ")
    )
    ifit = NA
    tryCatch({
      ifit = polr(
        formula(my_str)
        , data = data.frame(cbind(df1[[tab_in]], df2[[tab_out]]))
      )
      ifit$y = unlist(df2[[tab_out]][q_out])
      }
      , error=function(e){})
    ifit
  }
  , df1=scores
  , df2=data_tabs_f
)
surv2q$has_fit = sapply(
  surv2q$fit
  , function(ifit) class(ifit)=="polr"
)
# https://stats.stackexchange.com/questions/338904/measures-of-ordinal-classification-error-for-ordinal-regression
surv2q[surv2q$has_fit, "MSE"] = sapply( 
  surv2q[surv2q$has_fit, "fit"]
  , function(ifit) mean((fact2num(ifit$y) - fact2num(predict(ifit)))^2)
)
surv2q[surv2q$has_fit, "MSE_null"] = sapply(
  surv2q[surv2q$has_fit, "fit"]
  , function(ifit) {
    ipred = round(mean(fact2num(ifit$y)))
    mean((fact2num(ifit$y) - ipred)^2)
  })

surv2q[surv2q$has_fit, "MSE_adapted"] = sapply( 
  surv2q[surv2q$has_fit, "fit"]
  , function(ifit){
    # Balance responses
    y_counts = table(ifit$y)
    resids = fact2num(ifit$y) - fact2num(predict(ifit))
    weights = table(ifit$y)[paste0(fact2num(ifit$y))]
    (sum(resids^2*weights)/sum(weights))
  }
)
surv2q[surv2q$has_fit, "MSE_adapted_null"] = sapply(
  surv2q[surv2q$has_fit, "fit"]
  , function(ifit) {
    y_counts = table(ifit$y)
    weights =table(ifit$y)[paste0(fact2num(ifit$y))]
    ipred = round(mean(fact2num(ifit$y)))
    sum(weights*(fact2num(ifit$y) - ipred)^2) / sum(weights)
  })

surv2q[surv2q$has_fit,"dof"] = sapply(
  surv2q[surv2q$has_fit,"fit"] 
  , function(ifit) length(coef(ifit)) + length(ifit$zeta)
)
surv2q[surv2q$has_fit,"aic"] = sapply(surv2q[surv2q$has_fit,"fit"] , AIC)
surv2q$Question = as.vector(apply(
  surv2q, 1
  , function(irow){
    out_tab = paste0(irow["out_tab"])
    ilookup=lookups_all[[out_tab]]
    ilookup$.id = paste(out_tab, ilookup$ID, sep="_")
    out_q = paste0(irow["out_q"])
    istr = unlist(ilookup[which(ilookup$.id == out_q),"Question"])
    istr = paste0(istr, " ")
    istr
  }
))

# Out
surv2surv[,-3]
surv2score
surv2q[,-4]



###################
###### Output
###################



# Measures on THI
.surv2q_sub = surv2q[surv2q$out_tab=="THI",-c(4,5,7,9)]
.surv2q_join = inner_join(.surv2q_sub[.surv2q_sub$in_tab=="TFI",]
                          , .surv2q_sub[.surv2q_sub$in_tab=="TPFQ",]
                          , by="out_q", suffix=c("_TFI", "_TPFQ"))

# measure_tfi = .surv2q_join$aic_TFI; measure_tpfq = .surv2q_join$aic_TPFQ
# measure_tfi = .surv2q_join$MSE_adapted_TFI; measure_tpfq = .surv2q_join$MSE_adapted_TPFQ
measure_tfi = .surv2q_join$MSE_TFI; measure_tpfq = .surv2q_join$MSE_TPFQ
mean(measure_tfi < measure_tpfq) # How often is TFI better?
c("tfi"=mean(measure_tfi), "tpfq"=mean(measure_tpfq))
plot(measure_tfi, measure_tpfq
     , xlab="TFI MSE", ylab="TPFQ MSE"
     , main="MSE on Common Questions\nPredicting THI with TFI/TPFQ")
abline(0,1, col="darkred")


boxplot(.surv2score_sub$R2 ~ .surv2score_sub$in_tab)

boxplot(surv2q$MSE_adapted ~ surv2q$in_tab, ylim=c(0,1.5))


library(ggplot2)
.surv2score_sub = surv2score[(surv2score$in_tab != "THI") & (surv2score$out_tab != "THI"),]
.surv2score_sub$Prediction = .surv2score_sub$in_tab
levels(.surv2score_sub$Prediction) = c("TFI"="TFI to TPFQ", "THI"="THI", "TPFQ"="TPFQ to TFI")[levels(.surv2score_sub$Prediction)]
(ggplot(.surv2score_sub
       , aes(x=Prediction, y=R2)) 
  + geom_point()
  + labs(title="Subscales Prediction Accuracy" 
         , x ="Use Survey to Predict Survey"
         , y = "Explained Variance (%)")
  )



.surv2q_sub = surv2q # surv2q[(surv2q$out_tab == "TPFQ"),]
boxplot(.surv2q_sub$MSE ~ .surv2q_sub$in_tab,)

lookups_all$THI[c(8,19,24),]
