suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressPackageStartupMessages(library(IlluminaHumanMethylation450kmanifest))
suppressPackageStartupMessages(library(minfi))
suppressPackageStartupMessages(library(minfiData))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(CpGassoc))
suppressPackageStartupMessages(library(SmartSVA))
suppressPackageStartupMessages(library(isva))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(CpGFilter))
suppressPackageStartupMessages(library(structSSI))
suppressPackageStartupMessages(library(diptest))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(stringr))


Args <- commandArgs()
experiment <- Args[6]


#### load ref data (xReactiveProbes, ann450k) =====================
load("~/projects/ewas/analysis/reffile/ref.RData")
ann450k$refgene.pos <- sub(";.*","",ann450k$UCSC_RefGene_Group)
ann450k$refgene.pos <- ifelse(ann450k$refgene.pos=="",
                              "Non_gene",ann450k$refgene.pos)

#### set data directory and working directory =============================
outputDir <- "~/projects/ewas/analysis"
cmd <- paste0("mkdir -p ", outputDir, "/", experiment)
system(cmd)
baseDir <- "~/projects/ewas/data"
setwd(paste0(outputDir, "/", experiment))
dataset <- unlist(str_split(experiment, pattern = "_"))[1]
datadir <- paste0(baseDir, "/", dataset)


beta2m <- function(beta.values){
  log2(beta.values/(1 - beta.values))
}
m2beta <- function(m.values){
  2^m.values/(2^m.values + 1)
}

#### read in data and QC for samples ======================================
if(file.exists(paste0(datadir, "/SampleSheet.csv"))){
  targets <- read.metharray.sheet(datadir, pattern="SampleSheet.csv")
  rgSet <- read.metharray.exp(targets = targets)
  ## general information
  platform <- rgSet@annotation["array"]
  num_probes <- nrow(rgSet)
  col_name <- c()
  for(i in 1:ncol(rgSet)){
    samplename <- unlist(str_split(colnames(rgSet)[i], pattern = "_"))[1]
    col_name <- c(col_name, samplename)
  }
  colnames(rgSet) <- col_name
  ## read in phenotype information
  targets <- read.xlsx(paste0(baseDir, "/targets/", experiment, ".xlsx"))
  ## match phenotype info
  rgSet <- rgSet[, targets$Sample_Name]
  ## QC for samples
  detP <- detectionP(rgSet)
  keep <- colMeans(detP) < 0.01
  if(sum(keep) != nrow(targets)){
    cat("==========", sum(!keep)," sample was filted==========\n")
    rgSet <- rgSet[,keep]
    targets <- targets[keep,]
    detP <- detP[,keep]
  }
  sample_size <- ncol(rgSet)
  # Normalization
  mSetSq <- preprocessQuantile(rgSet)
  detP <- detP[rownames(mSetSq), ]
}else{
  filename1 <- paste0(datadir, "/", dataset, "_family.soft.gz")
  filename2 <- paste0(datadir, "/", dataset, "_family.soft")
  if(file.exists(filename1)){
    gset <- getGEO(filename = filename1)
  }else if(file.exists(filename2)){
    gset <- getGEO(filename = filename2)
  }else{
    stop("This dataset doesn't have soft file!")
  }
  ## general information
  platform <- gset@gpls[[1]]@header$title
  num_probes <- nrow(Table(gset@gsms[[1]]))
  ## beta values and m values
  bVals <- matrix(NA,nrow = nrow(Table(gset@gsms[[1]])),
                  ncol = length(gset@header$sample_id))
  rownames(bVals) <- Table(gset@gsms[[1]])$ID_REF
  colnames(bVals) <- gset@header$sample_id
  for (i in 1:length(gset@gsms)) {
    bVals[, i] <- Table(gset@gsms[[i]])[[2]]
  }
  if (all(is.na(bVals) | (bVals>=0 & bVals<=1))){
    mVals <- beta2m(bVals)
  } else {
    mVals <- bVals
    bVals <- m2beta(mVals)
  }
  
  ## P values (if exists)
  if (dim(Table(gset@gsms[[1]]))[2]>=3){
    detP <- matrix(NA,nrow = nrow(Table(gset@gsms[[1]])),
                   ncol = length(gset@header$sample_id))
    rownames(detP) <- rownames(bVals)
    colnames(detP) <- colnames(bVals)
    for (i in 1:length(gset@gsms)) {
      detP[, i] <- Table(gset@gsms[[i]])[[3]]
    }
  }
  
  ## read in phenotype information
  targets <- read.xlsx(paste0(baseDir, "/targets/", experiment, ".xlsx"))
  ## match phenotype info with values(beta, M, P) and probe info in ann450k
  bVals <- bVals[is.element(rownames(bVals), ann450k$Name), targets$Sample_Name]
  mVals <- mVals[is.element(rownames(mVals), ann450k$Name), targets$Sample_Name]
  if (exists("detP")){
    detP <- detP[is.element(rownames(detP), ann450k$Name), targets$Sample_Name]
    ## QC for samples
    keep <- colMeans(detP,na.rm=TRUE) < 0.01
    if(sum(keep) != nrow(targets)){
      cat("==========", sum(!keep)," sample was filtered==========\n")
      targets <- targets[keep,]
      detP <- detP[,keep]
      bVals <- bVals[, keep]
      mVals <- mVals[, keep]
    }
  }
  sample_size <- ncol(bVals)
  ## convert to GenomicRatioSet class
  mSetSq <- makeGenomicRatioSetFromMatrix(mat = bVals, pData = targets, what = "Beta")
  if(exists("detP")){detP <- detP[rownames(mSetSq), ]}
}

## summary
sample_info <- targets[, c("Sample_Name", "Sample_Group", "Rep_Design", "Cell_Type", "Tissue")]
write.csv(sample_info, file = paste0(experiment, "_sample_info.csv"), quote = FALSE, row.names = FALSE)
if(is.numeric(targets$Sample_Group)){
  index <- duplicated.data.frame(targets[, c("Cell_Type", "Tissue")])
  pheno_info <- targets[, c("Cell_Type", "Tissue")][!index, ]
  pheno_info$mean <- mean(targets$Sample_Group)
  pheno_info$size <- length(targets$Sample_Group)
}else{
  group_nums <- table(targets$Sample_Group)
  index <- duplicated.data.frame(targets[, c("Sample_Group", "Cell_Type", "Tissue")])
  pheno_info <- targets[, c("Sample_Group", "Cell_Type", "Tissue")][!index, ]
  pheno_info$size <- group_nums[pheno_info$Sample_Group]
}

#### QC for probes ========================================================
if (exists("detP")){
  keep <- rowSums(detP < 0.01,na.rm=TRUE) >= ncol(detP) * 0.5
  cat(paste(rep("=",10),collapse = ""),sum(keep==FALSE)," probes have bad quality ",
      paste(rep("=",10),collapse = ""),"\n")
  mSetSqFlt <- mSetSq[keep,]
} else {
  mSetSqFlt <- mSetSq
}

#filter out the probes from the X and Y chromosomes
keep <- !(featureNames(mSetSqFlt) %in%
            ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetSqFlt <- mSetSqFlt[keep,]
#filter out probes that have shown to be cross-reactive
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
mSetSqFlt <- mSetSqFlt[keep,]
#filter out the probes where common SNPs may affect the CpG or single base extension (SBE) site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG", "SBE"))
# the probes' number of experiment after filter
num_probes_filter <- nrow(mSetSqFlt)


#### get beta and M values =========================================
bVals <- getBeta(mSetSqFlt)

#Remove the lines with NA values
bVals <- na.omit(bVals)

#Adjust the beta values of 0/1
if(sum(bVals == 0) > 0){
  bVals[bVals == 0] <- min(bVals[bVals != 0])
}
if(sum(bVals == 1) > 0){
  bVals[bVals == 1] <- max(bVals[bVals != 1])
}

# M values
mVals <- beta2m(bVals)
## ICC analysis by cpgfilter
# only calculate icc when technical replicates exist
if(length(unique(targets$Rep_Design)) != length(targets$Rep_Design)){
  icc.m <- CpGFilterICC(mVals, targets$Rep_Design, logit.transform = FALSE)
  icc.b <- CpGFilterICC(bVals, targets$Rep_Design, logit.transform = FALSE)
}

#### use unique samples for the following analysis
index <- duplicated(targets$Rep_Design)
targets2 <- targets[!index, ]
mVals2 <- mVals[, targets2$Sample_Name]
bVals2 <- bVals[, targets2$Sample_Name]

#### SVA ==================================================================
mod <- model.matrix( ~ targets2$Sample_Group)

#for m value
Y.r <- t(resid(lm(t(mVals2) ~ targets2$Sample_Group)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
sva.res <- smartsva.cpp(mVals2, mod, mod0=NULL, n.sv=n.sv)

#### CpGassoc analysis ====================================================
if(nlevels(as.factor(targets2$Sample_Group)) == 2){
  cpgRes <- cpg.assoc(bVals2, as.numeric(as.factor(targets2$Sample_Group)), logit.transform=TRUE)
  cpgRes.sv <- cpg.assoc(bVals2, as.numeric(as.factor(targets2$Sample_Group)),
                         as.data.frame(sva.res$sv),logit.transform=TRUE)
}else{
  cpgRes <- cpg.assoc(bVals2, targets2$Sample_Group, logit.transform=TRUE)
  cpgRes.sv <- cpg.assoc(bVals2, targets2$Sample_Group,
                         as.data.frame(sva.res$sv),logit.transform=TRUE)
}
cpgResDf <- cpgRes$results
cpgResDf <- na.omit(cpgResDf)
write.csv(cpgResDf,file="cpgResDf.csv")

cpgResDf.sv <- cpgRes.sv$results
cpgResDf.sv <- na.omit(cpgResDf.sv)
write.csv(cpgResDf.sv,file="cpgResDf.sv.csv")


#### Present the datasets ==========================================
phenotype <- paste0(unique(targets$Sample_Group), collapse = ";")
celltype <- paste0(unique(targets$Cell_Type), collapse = ";")
tissue <- paste0(unique(targets$Tissue), collapse = ";")

## results before sva =================================
#Calculation of GIF, pi0, and significant count 
chisq <- qchisq(1-cpgResDf$P.value,1)
lambda <- median(chisq)/qchisq(0.5,1)
lambda <- round(lambda,3)
pi0 <- pi0est(cpgResDf$P.value,lambda = seq(0.05, 0.95, 0.05),pi0.method = "bootstrap")$pi0
pi0 <- round(pi0,3)

BH_sig_probes <- sum(p.adjust(cpgResDf$P.value, method = "fdr") < 0.05)
Holm_sig_probes <- sum(p.adjust(cpgResDf$P.value, method = "holm") < 0.05)
ST_sig_probes <- sum(qvalue(cpgResDf$P.value)$qvalues < 0.05)

## results after sva =================================
#Calculation of GIF, pi0, and significant count 
chisq.sv <- qchisq(1-cpgResDf.sv$P.value,1)
lambda.sv <- median(chisq.sv)/qchisq(0.5,1)
lambda.sv <- round(lambda.sv,3)
pi0.sv <- pi0est(cpgResDf.sv$P.value,lambda = seq(0.05, 0.95, 0.05),pi0.method = "bootstrap")$pi0
pi0.sv <- round(pi0.sv,3)

BH_pvalue <- p.adjust(cpgResDf.sv$P.value, method = "fdr")

ST_qvalue <- qvalue(cpgResDf.sv$P.value)

BH_sig_probes.sv <- sum( BH_pvalue < 0.05)
Holm_sig_probes.sv <- sum(p.adjust(cpgResDf.sv$P.value, method = "holm") < 0.05)
ST_sig_probes.sv <- sum(ST_qvalue$qvalues < 0.05)

data.qc <- data.frame(ID = experiment,
                      platform = platform,
                      phenotype = phenotype,
                      celltype = celltype,
                      tissue = tissue,
                      sample_size = sample_size,
                      num_probes = num_probes, 
                      num_probes_filter = num_probes_filter,
                      lambda_before_sva = lambda,
                      pi0_before_sva = pi0,
                      BH_sig_probes_before_sva = BH_sig_probes,
                      #ST_sig_probes_before_sva = ST_sig_probes,
                      Holm_sig_probes_before_sva = Holm_sig_probes,
                      lambda_after_sva = lambda.sv,
                      pi0_after_sva = pi0.sv,
                      BH_sig_probes_after_sva = BH_sig_probes.sv,
                      #ST_sig_probes_after_sva = ST_sig_probes.sv,
                      Holm_sig_probes_after_sva = Holm_sig_probes.sv
)
write.csv(data.qc,file="QC.csv", quote = FALSE, row.names = FALSE)

#### Covariates ==================================================================
covariate.name.list <- c("sd.b", "sd.m", "mean", "pos", "mad","dip",
                         "precision", "refgene.pos","cpg.loc", "chr", "dhs", "probe.type")
contin.cv <- c("sd.b", "sd.m", "mean", "mad", "dip", "precision", "pos")
cate.cv <- c("refgene.pos","cpg.loc","chr", "dhs", "probe.type")
if(length(unique(targets$Rep_Design)) != length(targets$Rep_Design)){
   covariate.name.list <- c(covariate.name.list, "icc.m", "icc.b")
   contin.cv <- c(contin.cv, "icc.m", "icc.b")
}

bVals2 <- bVals2[cpgResDf.sv$CPG.Labels, ]
mVals2 <- mVals2[cpgResDf.sv$CPG.Labels, ]
## Continuous covariates
sd.b <- apply(bVals2, 1, sd)
sd.m <- apply(mVals2, 1, sd)
mean <- apply(bVals2, 1, mean)

# mad: Median absolute deviation of beta values
mad <- apply(bVals2, 1, mad)

#dip: measure of unimodality in a sample
dip <- apply(bVals2, 1, dip)

#Precision: inverse precision parameter
precision <- 1/(mean*(1-mean)/sd.b^2-1)

pos <- ann450k$pos[match(cpgResDf.sv$CPG.Labels, rownames(ann450k))]

## Category covariates
refgene.pos <- ann450k$refgene_pos[match(cpgResDf.sv$CPG.Labels, rownames(ann450k))]
cpg.loc <- ann450k$Relation_to_Island[match(cpgResDf.sv$CPG.Labels, rownames(ann450k))]
chr <- ann450k$chr[match(cpgResDf.sv$CPG.Labels, rownames(ann450k))]
dhs <- ann450k$DHS[match(cpgResDf.sv$CPG.Labels, rownames(ann450k))]
dhs[!dhs == "TRUE"] <- "FALSE"   # replace "" in dhs vetor with "FALSE"
probe.type <- ann450k$Type[match(cpgResDf.sv$CPG.Labels, rownames(ann450k))]

# signs of the estimated effect size (for continuous and binary phenotype)
if(nlevels(as.factor(targets2$Sample_Group)) == 2 | is.numeric(targets2$Sample_Group)) {
  coef_df <- cpgRes.sv$coefficients[cpgResDf.sv$CPG.Labels, ]
  direction <- ifelse(coef_df$adj.intercept>0,"positive","negative")
  covariate.name.list <- c(covariate.name.list, "direction")
  cate.cv <- c(cate.cv, "direction")
}


##Combining covariates
cov.list <- list()
for(i in 1:length(covariate.name.list)){
  cov.list[[i]] <- get(covariate.name.list[i])
}
cov.df <- as.data.frame(cov.list)
colnames(cov.df) <- covariate.name.list
rownames(cov.df) <- rownames(bVals2)
write.csv(cov.df, file="Covariates.csv")

