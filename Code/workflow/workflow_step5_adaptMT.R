library(qvalue)
library(splines)
library(adaptMT)

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

### method AdaPT =============================
# degree of freedom
# res_glm$info
# res_glm$qvals

AdaPT_result <- list()
AdaPT_time <- list()
for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    X <- data.frame(x = cov.df[[i]])
    formulas <- "ns(x, df = 6)"
    starttime <- Sys.time()
    res_glm <- adapt_glm(x = X, pvals = cpgResDf.sv$P.value, pi_formulas = formulas, mu_formulas = formulas)
    AdaPT_result[[i]] <- res_glm$qvals
    endtime <- Sys.time()
    AdaPT_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
  }else if(is.element(i, cate.cv)){
    X <- data.frame(x = cov.df[[i]])
    formulas <- "x"
    starttime <- Sys.time()
    res_glm <- adapt_glm(x = X, pvals = cpgResDf.sv$P.value, pi_formulas = formulas, mu_formulas = formulas)
    AdaPT_result[[i]] <- res_glm$qvals
    endtime <- Sys.time()
    AdaPT_time[[i]] <- as.numeric(endtime - starttime, units = "secs")
  }
}

write.csv(AdaPT_result, file = "AdaPT_result.csv", row.names = FALSE, quote = FALSE)
write.csv(AdaPT_time, file = "AdaPT_time.csv", row.names = FALSE, quote = FALSE)

AdaPT_result <- as.data.frame(AdaPT_result)
AdaPT_time <- as.data.frame(AdaPT_time)


## compare probes' number ==================
AdaPT_probes <- list()
for(i in colnames(AdaPT_result)){
  AdaPT_probes[[i]] <- sum(AdaPT_result[[i]] < 0.05)
}
probes <- as.data.frame(AdaPT_probes)
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

write.csv(probes, file = "AdaPT_Discoveries.csv", quote = FALSE)
