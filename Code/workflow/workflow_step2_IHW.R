suppressPackageStartupMessages(library(IHW))
suppressPackageStartupMessages(library(qvalue))


Args <- commandArgs()
print(Args)
experiment <- Args[6]

outputDir <- "~/projects/ewas/analysis"
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
contin.cv <- c("sd.b", "sd.m", "mean", "mad","dip","precision","pos", "icc.b", "icc.m")
cate.cv <- c("refgene.pos","cpg.loc","chr", "dhs", "direction", "probe.type")
statistic.cv <- c("sd.b", "sd.m", "mean", "icc.b", "icc.m", "mad","dip","precision","direction")
CpGs.cv <- c("pos", "refgene.pos", "cpg.loc","chr", "dhs", "probe.type")

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

### method IHW ===============================
IHW_result <- list()
IHW_time <- list()
for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    covariate_type <- "ordinal"  ## for continous covariates
    starttime <- Sys.time()
    result <- ihw(cpgResDf.sv$P.value, cov.df[[i]], alpha = 0.05, covariate_type = covariate_type)
    endtime <- Sys.time()
  }else{
    covariate_type <- "nominal"  ## for category covariates
    starttime <- Sys.time()
    result <- ihw(cpgResDf.sv$P.value, as.factor(cov.df[[i]]), alpha = 0.05, covariate_type = covariate_type)
    endtime <- Sys.time()
  }
  IHW_result[[i]] <- adj_pvalues(result)
  IHW_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
}

write.csv(IHW_result, file = "IHW_result.csv", row.names = FALSE, quote = FALSE)
write.csv(IHW_time, file = "IHW_time.csv", row.names = FALSE, quote = FALSE)

IHW_result <- as.data.frame(IHW_result)
IHW_time <- as.data.frame(IHW_time)

## compare probes' number ==================
IHW_probes <- list()
for(i in colnames(IHW_result)){
  IHW_probes[[i]] <- sum(IHW_result[[i]] < 0.05)
}
probes <- as.data.frame(IHW_probes)
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

write.csv(probes, file = "IHW_Discoveries.csv", quote = FALSE)
