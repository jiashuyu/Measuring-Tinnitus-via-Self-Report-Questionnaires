library(ggplot2)
library(psych)
load("01B_out_efacfa_2021-03-03-1824-31.RData")

"THI analysis"
THI_results=results$THI

## EFA
KMO(data_tabs[['THI']])
cortest.bartlett(data_tabs[['THI']])
THI_EFA=THI_results$efa
THI_EFA3=THI_EFA[[3]]
# Unexplained Variance by Components
THI_EFA3$uniquenesses/(THI_EFA3$communalities+THI_EFA3$uniquenesses)
THI_EFA3$loadings
THI_EFA3$Vaccounted 
## rotation
### oblimin(output)
THI_EFA3$stat424_oblimin
### varimax
THI_EFA3$stat424_varimax
### target
THI_EFA3$stat424_target

## CFA
### full
legacy_cfa_full=THI_results$legacy_cfa$full
summary(legacy_cfa_full,standardized=TRUE,fit.measures = TRUE)
### reduced
legacy_cfa_reduced=THI_results$legacy_cfa$reduced
summary(legacy_cfa_reduced,standardized=TRUE,fit.measures = TRUE)


"TPFQ anaysis"
TPFQ_results=results$TPFQ

## EFA
KMO(data_tabs[['TPFQ']])
cortest.bartlett(data_tabs[['TPFQ']])
TPFQ_EFA=TPFQ_results$efa
TPFQ_EFA4=TPFQ_EFA[[4]]
# Unexplained Variance by Components
TPFQ_EFA4$uniquenesses/(TPFQ_EFA4$communalities+TPFQ_EFA4$uniquenesses)
TPFQ_EFA4$loadings
TPFQ_EFA4$Vaccounted 
## rotation
### oblimin
TPFQ_EFA4$stat424_oblimin
### varimax
TPFQ_EFA4$stat424_varimax
### target
TPFQ_EFA4$stat424_target

## CFA
### full
legacy_cfa_full=TPFQ_results$legacy_cfa$full
summary(legacy_cfa_full,standardized=TRUE,fit.measures = TRUE)
### reduced
legacy_cfa_reduced=TPFQ_results$legacy_cfa$reduced
summary(legacy_cfa_reduced,standardized=TRUE,fit.measures = TRUE)

## Delete Q3,5,8
data=read.csv('TPFQ.csv')
data=data[,-c(3,5,8)]
CorM=cor(data, use="complete.obs", method = "spearman")
model=fa(data,4, fm="ml",rotate="oblimin")
model$Vaccounted


"TFI anaysis"
TFI_results=results$TFI

## EFA
KMO(data_tabs[['TFI']])
cortest.bartlett(data_tabs[['TFI']])
TFI_EFA=TFI_results$efa
TFI_EFA8=TFI_EFA[[8]]
# Unexplained Variance by Components
TFI_EFA8$uniquenesses/(TFI_EFA8$communalities+TFI_EFA8$uniquenesses)
TFI_loading=TFI_EFA8$loadings
TFI_EFA8$Vaccounted 
## rotation
### oblimin
TFI_EFA8$stat424_oblimin
### varimax
TFI_EFA8$stat424_varimax
### target
TFI_EFA8$stat424_target
TFI_loading=TFI_EFA8$stat424_target$loadings

## CFA
### full
legacy_cfa_full=TFI_results$legacy_cfa$full
summary(legacy_cfa_full,standardized=TRUE,fit.measures = TRUE)
### reduced
legacy_cfa_reduced=TFI_results$legacy_cfa$reduced
summary(legacy_cfa_reduced,standardized=TRUE,fit.measures = TRUE)

## Delete Q3,5,8
data=read.csv('TFI.csv')
data=data[,-c(4,22)]
CorM=cor(data, use="complete.obs", method = "spearman")
model=fa(data,8, fm="ml",rotate="oblimin",maxit=2000)
model$Vaccounted

