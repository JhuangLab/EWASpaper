library(qvalue)
library(splines)

Args <- commandArgs()
print(Args)
experiment <- Args[6]


outputDir <- "/home/bailing/projects/ewas/analysis"
setwd(paste0(outputDir, "/", experiment))

cov.df <- read.csv("Covariates.csv", header = TRUE, row.names = 1)
covariate.name.list <- colnames(cov.df)
# replace NA with median value
for(i in 1:ncol(cov.df)){
  if(sum(is.na(cov.df[[i]])) != 0){
    cov.df[[i]][is.na(cov.df[[i]])] <- median(cov.df[[i]][!is.na(cov.df[[i]])])
  }
}

## classfication of covariates
contin.cv <- c("sd.b", "sd.m", "mean.b", "mad","dip","precision","pos", "icc.b", "icc.m")
cate.cv <- c("refgene.pos","cpg.loc","chr", "dhs", "direction", "probe.type")
statistic.cv <- c("sd.b", "sd.m", "var.b","var.m","mean.b", "icc.b", "icc.m",
                  "mad","dip","precision","direction")
CpGs.cv <- c("pos", "refgene.pos","cpg.loc","chr", "dhs", "probe.type")

for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    cov.df[[i]] <- as.numeric(cov.df[[i]])
  }else{
    cov.df[[i]] <- as.factor(cov.df[[i]])
  }
}
cpgResDf.sv <- read.csv("cpgResDf.sv.csv", header = TRUE, row.names = 1)
rownames(cpgResDf.sv) <- cpgResDf.sv$CPG.Labels
cpgResDf.sv <- cpgResDf.sv[rownames(cov.df), ]     

### method BH ===============================
starttime <- Sys.time()
BH_pvalue <- p.adjust(cpgResDf.sv$P.value, method = "fdr")
endtime <- Sys.time()
BH_time <- as.numeric(endtime - starttime, units = "secs")

### mthod ST ===============================
starttime <- Sys.time()
ST_qvalue <- qvalue(cpgResDf.sv$P.value)$qvalues
endtime <- Sys.time()
ST_time <- as.numeric(endtime - starttime, units = "secs")

### method CAMT ===============================
CAMT_result <- list()
CAMT_time <- list()
source("/home/bailing/projects/ewas/code/CAMT/camt.cor.func.R")
for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    X <- ns(cov.df[[i]], df = 6)
    starttime <- Sys.time()
    camt.obj <- camt(pvals = cpgResDf.sv$P.value, pi0.var = X, f1.var = X,
                     control.method = 'knockoff')
    endtime <- Sys.time()
    CAMT_result[[i]] <- camt.obj$fdr
    CAMT_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
  }else{
    starttime <- Sys.time()
    camt.obj <- camt(pvals = cpgResDf.sv$P.value, pi0.var = cov.df[[i]], f1.var = cov.df[[i]],
                     control.method = 'knockoff')
    endtime <- Sys.time()
    CAMT_result[[i]] <- camt.obj$fdr
    CAMT_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
  }
}

write.csv(CAMT_result, file = "CAMT_result.csv", row.names = FALSE, quote = FALSE)
write.csv(CAMT_time, file = "CAMT_time.csv", row.names = FALSE, quote = FALSE)

CAMT_result <- as.data.frame(CAMT_result)
CAMT_time <- as.data.frame(CAMT_time)

## compare probes' number ==================
CAMT_probes <- list()
for(i in colnames(CAMT_result)){
  CAMT_probes[[i]] <- sum(CAMT_result[[i]] < 0.05)
}
probes <- as.data.frame(CAMT_probes)
probes$BH <- sum(BH_pvalue < 0.05)
probes$ST <- sum(ST_qvalue < 0.05)

probes <- as.data.frame(t(probes))
colnames(probes) <- "discoveries"
probes$method <- rownames(probes)

group1 <- c("BH", "ST")
group2 <- rownames(probes)[rownames(probes) %in% statistic.cv]
group3 <- rownames(probes)[rownames(probes) %in% CpGs.cv]
method_rank <- c(group1, group2, group3)
probes <- probes[method_rank, ]
# probes <- probes[order(probes$discoveries), ]

write.csv(probes, file = "CAMT_Discoveries.csv", quote = FALSE)
