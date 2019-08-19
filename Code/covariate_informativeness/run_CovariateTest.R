Args <- commandArgs()
print(Args)
experiment <- Args[6]

setwd(paste0("~/projects/ewas/analysis/", experiment))
options(stringsAsFactors = FALSE)
# build a new directory to save results
cmd <- paste0("mkdir -p CovariatesTest")
system(cmd)

pvalue_df <- read.csv("cpgResDf.sv.csv", header = TRUE, row.names = 1)
covariates_df <- read.csv("Covariates.csv", header = TRUE, row.names = 1)
identical(pvalue_df$CPG.Labels, rownames(covariates_df))
covariates_df <- covariates_df[pvalue_df$CPG.Labels, ]
identical(pvalue_df$CPG.Labels, rownames(covariates_df))

# make sure the type of covariates is proporriate
covariate.name.list <- colnames(covariates_df)
# replace NA with median value
for(i in 1:ncol(covariates_df)){
  if(sum(is.na(covariates_df[[i]])) != 0){
    covariates_df[[i]][is.na(covariates_df[[i]])] <- median(covariates_df[[i]][!is.na(covariates_df[[i]])])
  }
}
## classfication of covariates
contin.cv <- c("sd.b", "sd.m", "mean", "mad","dip","precision","pos", "icc.b", "icc.m")
cate.cv <- c("refgene.pos","cpg.loc","chr", "dhs", "direction", "probe.type")
for(i in covariate.name.list){
  if(is.element(i, contin.cv)){
    covariates_df[[i]] <- as.numeric(covariates_df[[i]])
  }else{
    covariates_df[[i]] <- as.factor(covariates_df[[i]])
  }
}

# rank CpGs by chromosome position
library(tidyverse)
covariates <- tbl_df(covariates_df)
covariates$CPG.Label <- rownames(covariates_df)
covariates$chr <- factor(covariates$chr, levels = paste0("chr", 1:22))
ranked_CPGs <- arrange(covariates, chr, pos)$CPG.Label
pvalue_df <- pvalue_df[match(ranked_CPGs, pvalue_df$CPG.Labels), ]
covariates_df <- covariates_df[ranked_CPGs, ]
# find min auto correlation
x <- acf(pvalue_df$P.value, type = "correlation")
by_n <- x$lag[which.min(abs(x$acf))]
sprintf("downsample was carried out by every %.0f CpGs!", by_n)

# downsample
index <- seq.int(1, nrow(pvalue_df), by = by_n)
pvalue_df <- pvalue_df[index, ]
p_value <- pvalue_df$P.value
covariates_df <- covariates_df[index, ]

source("~/projects/ewas/code/CovariateTest/CovariateTest.R")
stat.o.a <- list()
stat.o <- list()
stat.p <- list()
p.value <- list()
stat.o1 <- list()
stat.p1 <- list()
p.value1 <- list()
stat.o2 <- list()
stat.p2 <- list()
p.value2 <- list()
x.cut.optim <- list()
p.cut.optim <- list()

for(i in covariate.name.list){
  print(paste0("start to deal with covariate ", i))
  if(is.numeric(covariates_df[[i]])){
    if(! i %in% c("icc.m", "icc.b")){
      output <- CovTest.n(pvalue = p_value, covariate = covariates_df[[i]])
    }else{
      # add small disturbance to zero values
      icc <- covariates_df[[i]]
      icc[icc == 0] <- runif(sum(icc == 0), 0, min(icc[icc != 0]))
      output <- CovTest.n(pvalue = p_value, covariate = icc)
    }
  }else if(is.factor(covariates_df[[i]])){
    output <- CovTest.c(pvalue = p_value, covariate = covariates_df[[i]])
  }
  stat.o.a[[i]] <- output$stat.o.a
  stat.o[[i]] <- output$stat.o
  stat.p[[i]] <- output$stat.p
  p.value[[i]] <- output$p.value
  stat.o1[[i]] <- output$stat.o1
  stat.p1[[i]] <- output$stat.p1
  p.value1[[i]] <- output$p.value1
  stat.o2[[i]] <- output$stat.o2
  stat.p2[[i]] <- output$stat.p2
  p.value2[[i]] <- output$p.value2
  x.cut.optim[[i]] <- output$x.cut.optim
  p.cut.optim[[i]] <- output$p.cut.optim
  print(paste0(i, " complete !"))
}

stat.o.a <- as.data.frame(stat.o.a)
stat.p <- as.data.frame(stat.p)
p.cut.optim <- data.frame(measure = rep(names(p.cut.optim), unlist(lapply(p.cut.optim, length))),
                          cutoff = unlist(lapply(p.cut.optim, names), use.names = FALSE),
                          p.cut.optim = unlist(p.cut.optim, use.names = FALSE))

# output stat.o, p.value, x.cut.optim and stat.o1, ect,into one dataframe
result <- data.frame(stat.o = unlist(stat.o),
                     p.value = unlist(p.value))

result_numeric <- data.frame(stat.o1 = unlist(stat.o1),
                             p.value1 = unlist(p.value1),
                             stat.o2 = unlist(stat.o2),
                             p.value2 = unlist(p.value2),
                             x.cut.optim = unlist(x.cut.optim))
final <- merge(result, result_numeric, by.x = 0, by.y = 0, all = TRUE)
final$dataset <- experiment

write.csv(final, file = "CovariatesTest/CovariatesTest_result.csv", quote = FALSE, row.names = FALSE)
write.csv(stat.o.a, file = "CovariatesTest/stat.o.a.csv", quote = FALSE, row.names = FALSE)
write.csv(stat.p, file = "CovariatesTest/stat.p.csv", quote = FALSE, row.names = FALSE)
write.csv(p.cut.optim, file = "CovariatesTest/p.cut.optim.csv", quote = FALSE, row.names = FALSE)