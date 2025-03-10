library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(openxlsx)

# xReactiveProbes object
xReactiveProbesDir <- "~/projects/ewas/data/ref/48639-non-specific-probes-Illumina450k.xlsx"
xReactiveProbes <- read.xlsx(xReactiveProbesDir,sheet = "nonspecific cg probes")


# ann450k object
ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

## save objects =====================
save(xReactiveProbes, ann450k, file = "~/projects/ewas/analysis/reffile/ref.RData")
