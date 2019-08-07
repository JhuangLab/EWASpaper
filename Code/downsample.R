######### define the true positive probes =====================================
experiment <- "GSE80261"
setwd(paste0("/home/bailing/projects/ewas/analysis/downsample/", experiment))
options(stringsAsFactors = FALSE)
data_file <- paste0("/home/bailing/projects/ewas/analysis/", experiment, "/cpgResDf.sv.csv")
raw_pvalue <- read.csv(data_file, header = TRUE, row.names = 1)
true_positive <- raw_pvalue$CPG.Labels[p.adjust(raw_pvalue$P.value, method = "bonferroni") < 0.05]
write.table(true_positive, file = "defined_true_positive_probes.txt", quote = FALSE, row.names = FALSE)

######### resampling at various sample sizes===================================
# run Rscript command in bash, for example "Rscript downsample.R GSE80261 10"

## setup sample size
Args <- commandArgs()
print(Args)
experiment <- Args[6]
sample_size <- as.numeric(Args[7])

## load beta values and phenotype information
load(paste0("/home/bailing/projects/ewas/analysis/", experiment, "/bvals.rdata"))
library(readxl)
targets <- read_xlsx(paste0("/home/bailing/projects/ewas/data/targets/", experiment, ".xlsx"))
if(is.numeric(targets$Sample_Group)){
  print("This dataset has continous category !")
}else{
  groups <- unique(targets$Sample_Group)
  group1_index <- which(targets$Sample_Group == groups[1])
  group2_index <- which(targets$Sample_Group == groups[2])
}

## read in the defined true positive probes list
true_positive <- read.table(paste0("/home/bailing/projects/ewas/analysis/downsample/", experiment, "/defined_true_positive_probes.txt"),
                            header = TRUE)
true_positive <- true_positive$x

## convert b value to m values
beta2m <- function(beta.values){
  log2(beta.values/(1 - beta.values))
}

## load 450k probes annoation information
load("/home/bailing/projects/ewas/analysis/reffile/ref.RData")
ann450k$refgene_pos <- sub(";.*","",ann450k$UCSC_RefGene_Group)
ann450k$refgene_pos <- ifelse(ann450k$refgene_pos=="",
                              "Non_gene",ann450k$refgene_pos)

## setup function to calculate Power
cal_Power <- function(x){
  Power <- length(intersect(x, true_positive)) / length(true_positive)
  return(Power)
}

## setup empty data frame to save results
cols_name <- c("BH", "ST", "sd.b", "sd.m", "mean.b", "pos", "mad",
               "dip", "precision", "refgene.pos","cpg.loc", "chr", "dhs", 
               "probe.type", "direction")
sensitivity <- matrix(NA, nrow = 100, ncol = 17, dimnames = list(NULL, cols_name))


library(minfi)
library(qvalue)
library(CpGassoc)
library(SmartSVA)
library(diptest)
library(splines)
# take method CAMT as an example
source("/home/bailing/projects/ewas/code/CAMT/camt.cor.func.R")

## sampling and calculation
for(j in 1:100){
  print(paste0("Replication ", j, " ============================================="))
  tryCatch({
    # randomly extract size i samples to set up a subset
    if(is.numeric(targets$Sample_Group)){
      targets_subset <- targets[sample(1:nrow(targets), sample_size), ]
    }else{
      targets_subset <- targets[c(sample(group1_index, sample_size), sample(group2_index, sample_size)), ]
    }
    bVals_matrix <- bVals[, targets_subset$Sample_Name]
    
    # convert to GenomicRatioSet class
    mSetSq <- makeGenomicRatioSetFromMatrix(mat = bVals_matrix, pData = targets_subset, what = "Beta")
    
    ## QC for probes ------------------------------------------------
    if (exists("detP")){
      detP_matrix <- detP[rownames(mSetSq), targets_subset$Sample_Name]
      keep <- rowSums(detP_matrix < 0.01,na.rm=TRUE) >= ncol(detP_matrix) * 0.5
      cat(paste(rep("=",10),collapse = ""),sum(keep==FALSE)," probes have bad quality ",
          paste(rep("=",10),collapse = ""),"\n")
      mSetSqFlt <- mSetSq[keep,]
    } else {
      mSetSqFlt <- mSetSq
    }
    # filter out the probes from the X and Y chromosomes
    keep <- !(featureNames(mSetSqFlt) %in%
                ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
    mSetSqFlt <- mSetSqFlt[keep,]
    # filter out probes that have shown to be cross-reactive
    keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
    mSetSqFlt <- mSetSqFlt[keep,]
    # filter out the probes where common SNPs may affect the CpG or single base extension (SBE) site
    mSetSqFlt <- dropLociWithSnps(mSetSqFlt, snps = c("CpG", "SBE"))
    
    ## get beta values and M values ---------------------------------
    bVals_matrix <- getBeta(mSetSqFlt)
    bVals_matrix <- na.omit(bVals_matrix)
    #Adjust the beta values of 0/1
    if(sum(bVals_matrix == 0) > 0){
      bVals_matrix[bVals_matrix == 0] <- min(bVals_matrix[bVals_matrix != 0])
    }
    if(sum(bVals_matrix == 1) > 0){
      bVals_matrix[bVals_matrix == 1] <- max(bVals_matrix[bVals_matrix != 1])
    }
    
    mVals_matrix <- beta2m(bVals_matrix)
    
    ## SVA based on M values ----------------------------------------
    mod <- model.matrix( ~ targets_subset$Sample_Group)
    Y.r <- t(resid(lm(t(mVals_matrix) ~ targets_subset$Sample_Group)))
    n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
    sva.res <- smartsva.cpp(mVals_matrix, mod, mod0=NULL, n.sv=n.sv)
    
    ## calculate p values -------------------------------------------
    if(is.numeric(targets_subset$Sample_Group)){
      cpgRes.sv <- cpg.assoc(bVals_matrix, targets_subset$Sample_Group,
                             as.data.frame(sva.res$sv),logit.transform=TRUE)
    }else{
      cpgRes.sv <- cpg.assoc(bVals_matrix, as.numeric(as.factor(targets_subset$Sample_Group)),
                             as.data.frame(sva.res$sv),logit.transform=TRUE)
    }
    cpgResDf.sv <- na.omit(cpgRes.sv$results)
    bVals_matrix <- bVals_matrix[cpgResDf.sv$CPG.Labels, ]
    mVals_matrix <- mVals_matrix[cpgResDf.sv$CPG.Labels, ]
    
    #### obtain different covariates --------------------------------
    
    ## obtain the 450k probes information as covariates -------------
    # chromosome position, like 53468112,37459206,etc.
    pos <- ann450k$pos[match(rownames(bVals_matrix), rownames(ann450k))]
    # refgene position, like Body,3'UTR,etc.
    refgene.pos <- ann450k$refgene_pos[match(rownames(bVals_matrix), rownames(ann450k))]
    # CpG location, like N_Shore,OpenSea,etc.
    cpg.loc <- ann450k$Relation_to_Island[match(rownames(bVals_matrix), rownames(ann450k))]
    # chromosome number, like chr16,chr1,etc.
    chr <- ann450k$chr[match(rownames(bVals_matrix), rownames(ann450k))]
    # DNAase hypersensitive site
    dhs <- ann450k$DHS[match(rownames(bVals_matrix), rownames(ann450k))]
    dhs[!dhs == "TRUE"] <- "FALSE"   # replace "" in DHS vetor with "FALSE"
    # illumina 450k array probe type, two types in total.
    probe_type <- ann450k$Type[match(rownames(bVals_matrix), rownames(ann450k))]
    
    ## obtain statistic variable as covariates ----------------------
    sd.b <- apply(bVals_matrix, 1, sd)
    sd.m <- apply(mVals_matrix, 1, sd)
    mean.b <- apply(bVals_matrix, 1, mean)
    mad <- apply(bVals_matrix, 1, mad)
    dip <- apply(bVals_matrix,1,dip)
    precision <- 1/(mean.b*(1-mean.b)/sd.b^2-1)
    coef_df <- cpgRes.sv$coefficients[cpgResDf.sv$CPG.Labels, ]
    direction <- ifelse(coef_df$adj.intercept>0,"positive","negative")
    
    #### different methods to adjust p values -----------------------
    # method BH
    probes <- cpgResDf.sv$CPG.Labels[p.adjust(cpgResDf.sv$P.value, method = "fdr") < 0.05]
    sensitivity[j, "BH"] <- cal_Power(probes)
    # method ST
    probes <- cpgResDf.sv$CPG.Labels[qvalue(cpgResDf.sv$P.value)$qvalues < 0.05]
    sensitivity[j, "ST"] <- cal_Power(probes)
    # continuous variates
    for(s in c("sd.b", "sd.m", "mean.b", "mad", "dip", "precision", "pos")){
      camt.obj <- camt(pvals = cpgResDf.sv$P.value, pi0.var = ns(get(s), df = 6), f1.var = ns(get(s), df = 6),
                       control.method = 'knockoff')
      probes <- cpgResDf.sv$CPG.Labels[camt.obj$fdr < 0.05]
      sensitivity[j, s] <- cal_Power(probes)
    }
    # category variates
    for(k in c("direction", "refgene.pos", "cpg.loc", "chr", "dhs", "probe.type")){
      camt.obj <- camt(pvals = cpgResDf.sv$P.value, pi0.var = as.factor(get(k)), f1.var = as.factor(get(k)),
                       control.method = 'knockoff')
      probes <- cpgResDf.sv$CPG.Labels[camt.obj$fdr < 0.05]
      sensitivity[j, k] <- cal_Power(probes)
    }
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

write.csv(sensitivity, file = paste0("/home/bailing/projects/ewas/analysis/downsample/", experiment, "/CAMT/power_", sample_size, ".csv"),
          row.names = TRUE, quote = FALSE)

