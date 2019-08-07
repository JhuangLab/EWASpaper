library(ggplot2)
myqqplot <- function(experiment){
  ## tidy p value data
  fn <- paste0("/home/bailing/projects/ewas/analysis/", experiment, "/")
  cpgResDf <- read.csv(paste0(fn, "cpgResDf.csv"), header = TRUE, row.names = 1)
  cpgResDf$label <- "Before SVA"
  cpgResDf$P.value[cpgResDf$P.value == 0] <- min(cpgResDf$P.value[cpgResDf$P.value != 0])
  cpgResDf.sv <- read.csv(paste0(fn, "cpgResDf.sv.csv"), header = TRUE, row.names = 1)
  cpgResDf.sv$label <- "After SVA"
  cpgResDf.sv$P.value[cpgResDf.sv$P.value == 0] <- min(cpgResDf.sv$P.value[cpgResDf.sv$P.value != 0])

  ## calculate the expected p value
  cpgResDf$exp.pvalue <- (rank(cpgResDf$P.value, ties.method="first")-.5)/(length(cpgResDf$P.value)+1)
  cpgResDf.sv$exp.pvalue <- (rank(cpgResDf.sv$P.value, ties.method="first")-.5)/(length(cpgResDf.sv$P.value)+1)
  
  ## merge data ready for plot
  df <- rbind(cpgResDf[, c("P.value", "label", "exp.pvalue")], cpgResDf.sv[, c("P.value", "label", "exp.pvalue")])
  df_plot <- unique(data.frame(P.value = round(-log10(df$P.value), 2),
                          label = df$label,
                          exp.pvalue = round(-log10(df$exp.pvalue), 2)))
  ## ggplot
  ggplot(data = df_plot, aes(x = exp.pvalue, y = P.value)) + geom_point(aes(color = label)) + geom_abline(intercept = 0, slope = 1) + xlab("Expected.-log10(p-value)") + ylab("Observed.-log10(p-value)") + theme_classic()
}

library(readxl)
id_conversion <- as.data.frame(read_xlsx("/home/bailing/projects/ewas/doc/id_conversion.xlsx", sheet = "Sheet1"))
dataset <- id_conversion$Raw_ID

library(futile.logger)
myplot <- list()
for(s in dataset){
  flog.info("start plotting for %s!", s)
  myplot[[s]] <- myqqplot(s)
  flog.info("done for %s!", s)
}
names(myplot) <- id_conversion$New_ID[match(names(myplot), id_conversion$Raw_ID)]

save(myplot, file = "/home/bailing/projects/ewas/analysis/qqplot.RData")
