library(qvalue)
library(splines)
library(FDRreg)

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
statistic.cv <- c("sd.b", "sd.m", "mean.b", "icc.b", "icc.m",
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

### method FDRreg ============================
## transform pvalues to z-score
ztransform <- function(p.value, tol = 1E-15) {
  # Transform p-values to z-scores
  p.value[p.value <= tol] <- tol
  p.value[p.value >= 1 - tol] <- 1 - tol
  z <- -qnorm(p.value)
  return(z)
}
zscore <- ztransform(cpgResDf.sv$P.value)

FDRreg_result <- list()
FDRreg_time <- list()
for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    X <- ns(cov.df[[i]], df = 6)
    starttime <- Sys.time()
    result <- FDRreg(zscore, X, nulltype = 'theoretical')
    FDRreg_result[[i]] <- result$FDR
    endtime <- Sys.time()
    FDRreg_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
    print(paste0(i, " done !"))
  }else{
    X <- model.matrix( ~ cov.df[[i]])[, -1]
    starttime <- Sys.time()
    result <- FDRreg(zscore, as.matrix(X), nulltype = 'theoretical')
    FDRreg_result[[i]] <- result$FDR
    endtime <- Sys.time()
    FDRreg_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
    print(paste0(i, " done !"))
  }
}

write.csv(FDRreg_result, file = "FDRreg_result.csv", row.names = FALSE, quote = FALSE)
write.csv(FDRreg_time, file = "FDRreg_time.csv", row.names = FALSE, quote = FALSE)

FDRreg_result <- as.data.frame(FDRreg_result)
FDRreg_time <- as.data.frame(FDRreg_time)


## compare probes' number ==================
FDRreg_probes <- list()
for(i in colnames(FDRreg_result)){
  FDRreg_probes[[i]] <- sum(FDRreg_result[[i]] < 0.05)
}
probes <- as.data.frame(FDRreg_probes)
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

write.csv(probes, file = "FDRreg_Discoveries.csv", quote = FALSE)
