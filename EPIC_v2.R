# DNA Methylation Pipeline for Illumina HumanMethylationEPIC BeadChip Data - NOT IDAT files but BeadStudio output
# WHEN NO. OF SAMPLES < 411 FOR EPIC DATA
# Developed by Suparna Kumar
# For the analysis of RA vs. Healthy whole (peripheral) blood samples

#### CITATIONS: #####

## For package "minfi"
# Martin J. Aryee, Andrew E. Jaffe, Hector Corrada-Bravo, Christine Ladd-Acosta,
# Andrew P. Feinberg, Kasper D. Hansen, and Rafael A. Irizarry, Minfi: A flexible
# and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation
# microarrays, Bioinformatics 30 (2014), no. 10, 1363–1369.

## For package "missMethyl"
# Jovana Maksimovic, Lavinia Gordon, and Alicia Oshlack, Swan: Subset-quantile within
# array normalization for illumina infinium humanmethylation450 beadchips, Genome
# Biology 13 (2012), no. 6, R44.

# Belinda Phipson and Alicia Oshlack, Diffvar: A new method for detecting differential
# variability with application to methylation in cancer and aging, Genome Biology 15 (2014), no. 9, 465.

# Jovana Maksimovic, Johann A Gagnon-Bartsch, Terence P Speed, and Alicia Oshlack,
# Removing unwanted variation in a differential methylation analysis of illumina
# humanmethylation450 array data, Nucleic acids research (2015), gkv526.

# Belinda Phipson and Jovana Maksimovic and Alicia Oshlack, missMethyl: an R
# package for analysing methylation data from Illuminas HumanMethylation450 platform, Bioinformatics
# (2015), btv560.

## For package "IlluminaHumanMethylation450kmanifest"
# Kasper Daniel Hansen and Martin Aryee (2012). IlluminaHumanMethylation450kmanifest: Annotation for
# Illumina's 450k methylation arrays. R package version 0.4.0.

## For package "IlluminaHumanMethylation450kanno.ilmn12.hg19"
# Kasper Daniel Hansen (2016). IlluminaHumanMethylation450kanno.ilmn12.hg19: Annotation for Illumina's
# 450k methylation arrays. R package version 0.6.0.

## For package "limma"
# Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers
# differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research
# 43(7), e47.

# Phipson, B, Lee, S, Majewski, IJ, Alexander, WS, and Smyth, GK (2016). Robust
# hyperparameter estimation protects against hypervariable genes and improves power
# to detect differential expression. Annals of Applied Statistics 10(2), 946–963.




################  START PIPELINE ####################


rm(list=ls())

# Loading all required libraries
library(missMethyl)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

## Function to obtain the raw p value and t statistics from the linear model fit produced by "varFit"
get_p_table <- function(fit, coef, number) {
  out <- data.frame(AvgVar = fit$AvgVar, LogVarRatio = fit$LogVarRatio[, coef], DiffLevene = fit$coeff[, coef], t = fit$t[, coef], P.Value = fit$p.value[, coef])
  out[1:number, ]
}

## Function to perform QC

performQC <- function(mSet, filename){
  # Quality control plot uses the log median intensity in both the methylated (M)
  # and unmethylated (U) channels. When plotting these two medians against each other,
  # it has been observed that good samples cluster together, while failed samples
  # tend to separate and have lower median intensities.
  qc <- addQC(mSet, getQC(mSet)) # adding phenotype data
  
  minfi_meth <- getMeth(mSet) # get methylated intensities
  minfi_unmeth <- getUnmeth(mSet) # get unmethylated intensities
  # Shorten column names to include only sample names
  colnames(minfi_meth) <- sub("\\_.*", "", colnames(minfi_meth))
  colnames(minfi_unmeth) <- sub("\\_.*", "", colnames(minfi_unmeth))
  log_minfi_meth <- log(minfi_meth)
  log_minfi_unmeth <- log(minfi_unmeth)
  
  #### Save in PDF ####
  pdf(file = filename, width = 10, height = 10)
  boxplot(log_minfi_meth, las = 2, cex.axis = 0.8, main = "Methylated")
  boxplot(log_minfi_unmeth, las = 2, cex.axis = 0.8, main = "Unmethylated")
  
#  plotQC(qc)
  badSampleCutoff = 10.5
  meds <- (qc$mMed + qc$uMed)/2 # mean of M & U medians per sample
  whichBad <- which((meds < badSampleCutoff))
  plot(qc$mMed, qc$uMed, xlim = c(8, 14), ylim = c(8, 14), 
       xaxt = "n", yaxt = "n", xlab = "Meth median intensity (log2)", 
       ylab = "Unmeth median intensity (log2)", col = ifelse(1:nrow(qc) %in% 
                                                               whichBad, "red", "black"))
  axis(side = 1, at = c(9, 11, 13))
  axis(side = 2, at = c(9, 11, 13))
  abline(badSampleCutoff * 2, -1, lty = 2)
  if (length(whichBad) > 0) {
    text(qc$mMed[whichBad], qc$uMed[whichBad] - 0.25, labels = whichBad, 
         col = "red")
  }
  legend("topleft", legend = c("good", "bad, with sample index"), 
         pch = 1, col = c("black", "red"), bty = "n")
  dev.off()
  
  badSampleNames <- mSet@colData@rownames[whichBad]
  sink(file = "/Users/sups/Downloads/R_Prog/COV/bad_samples.txt")
  print(badSampleNames)
  sink()
}

#### Step 1: Read the Illumina normalised EPIC intensities #####

# Import the EPIC data (methylated, unmethylated, detection p-values) into a table
intensities <- read.table("/Users/sups/Downloads/R_Prog/toy2.tsv", header = TRUE, sep = "\t")

# Create a table with only detection p-values to carry out sample filtering
include_indices <- grep("*_Pval", colnames(intensities))
det_table <- intensities[, include_indices]

#### Step 2: Filter out Y chromosome signals for females #####
# Only select those intensities that are lesser than the cut-off
# Detection p-value represents the confidence measure for the beta value of the sample
keep <- rowSums(det_table) < 0.01 # where 0.01 is the selected p-value
intensities <- intensities[keep,]
rm(list=c("det_table", "keep")) # Remove unnecessary variables

#### Step 3: Convert the M & U intensities into a MethylSet #####

# Converting to MethylSet object
manifest <- IlluminaHumanMethylationEPICmanifest
manifestName <- "IlluminaHumanMethylationEPICmanifest"

um_index <- grep("*_Un*", colnames(intensities))
unmethData <- intensities[, um_index]
rownames(unmethData) <- intensities$ID_REF
unmethData <- data.matrix(unmethData)
m_index <- grep("*_Methylated_signal*", colnames(intensities))
methData <- intensities[, m_index]
rownames(methData) <- intensities$ID_REF
methData <- data.matrix(methData)

colnames(methData)<- sub("\\_.*", "", colnames(methData))
colnames(unmethData)<- sub("\\_.*", "", colnames(unmethData))
intensities <- intensities[, -include_indices]
rownames(intensities) <- intensities$ID_REF
intensities <- intensities[, -1]
intensities <- SummarizedExperiment(assays=SimpleList("Meth" = methData, "Unmeth" = unmethData))

mSet <- MethylSet(Meth = methData, Unmeth = unmethData, colData = colData(intensities), annotation = c("array" = "IlluminaHumanMethylationEPIC", "annotation" = "ilm10b4.hg19"))
featureNames(mSet) <- intensities@NAMES
sampleNames(mSet) <- intensities@colData@rownames

mSet@preprocessMethod[["rg.norm"]] = "Illumina"
mSet@preprocessMethod[["minfi"]] = as.character(packageVersion("minfi"))
mSet@preprocessMethod[["manifest"]] = as.character(packageVersion(manifestName))
mSet@preprocessMethod = mSet@preprocessMethod[-1]

rm(manifest, manifestName, include_indices, m_index, um_index, methData, unmethData) # Remove unnecessary variables

# QC of the mSet
performQC(mSet, filename = "/Users/sups/Downloads/R_Prog/COV/Boxplot_Intensities_Outliers_Raw.pdf")

#### Step 4: SWAN normalization ####
mSetSw <- SWAN(mSet,verbose=FALSE)
# QC of SWAN normalized mSet
performQC(mSetSw, filename = "/Users/sups/Downloads/R_Prog/COV/Boxplot_Intensities_Outliers_SWAN.pdf")

# Plotting density distribution of beta values before and after using SWAN.
pdf(file = "/Users/sups/Downloads/R_Prog/COV/SWAN_Density.pdf", width = 10, height = 10)
par(mfrow=c(1,1), cex=0.50)
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")
dev.off()

# Extract beta and M-values from the SWAN normalised data.
# We prefer to add an offset to the methylated and unmethylated intensities
# when calculating M-values, hence we extract the methylated and unmethylated channels
# separately and perform our own calculation.
meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSw)

colnames(meth) <- sub("\\_.*", "", colnames(meth))
colnames(unmeth) <- sub("\\_.*", "", colnames(unmeth))
meth_colnames <- colnames(meth)
unmeth_colnames <- colnames(unmeth)
# Save the methylated sample names in a text file
sink(file = "/Users/sups/Downloads/R_Prog/COV/meth_names.txt")
print(meth_colnames)
sink()

# Save the unmethylated sample names in a text file
sink(file = "/Users/sups/Downloads/R_Prog/COV/unmeth_names.txt")
print(unmeth_colnames)
sink()

# Plot MDS (multi-dimensional scaling) of RA and normal samples.
# This is a good check to make sure samples cluster together according to their type.
pdf(file = "/Users/sups/Downloads/R_Prog/COV/Clustering.pdf", width = 10, height = 10)
par(mfrow=c(1,1))
plotMDS(Mval, labels=targets$Sample_Name, col=as.integer(factor(targets$Sample_Group)))
legend("bottom",legend=c("Healthy","RA"),pch=16,cex=1.2,col=1:2)

tree <- hclust(dist(t(Mval)))
plot(tree)

dev.off()
# We test for differential methylation using the *limma* package which
# employs an empirical Bayes framework based on Guassian model theory.
# We need to set up the design matrix. The most straightforward is directly
# from the targets file.

# Create experimental design using limma
design <- model.matrix(~targets$Sample_Group)

# Now we can test for differential methylation using the lmFit and eBayes
# functions from limma. As input data we use the matrix of M-values.
fit.reduced <- lmFit(Mval,design)
fit.reduced <- eBayes(fit.reduced)

summary(decideTests(fit.reduced))
top<-topTable(fit.reduced,coef=2) # coef is the last column index in the decideTests table

cpgs <- rownames(top)
pdf(file = "/Users/sups/Downloads/R_Prog/COV/4DMR.pdf", width = 10, height = 10)
par(mfrow=c(2,2))
for(i in 1:4){
  stripchart(beta[rownames(beta)==cpgs[i],]~design[,2],method="jitter",
             group.names=c("Healthy","RA"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
             vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
  title(cpgs[i],cex.main=1.5)
}
dev.off()

#### Step 5: Batch normalization ####

# 450k array studies are subject to unwanted technical variation such as batch
# effects and other, often unknown, sources of variation. The adverse effects of
# unwanted variation have been extensively documented in gene expression array
# studies and have been shown to be able to both reduce power to detect true differences
# and to increase the number of false discoveries.

# Negative control features are probes/genes/etc. that are known a priori to not truly
# be associated with the biological factor of interest, but are affected by unwanted
# variation. For example, in a microarray gene expression study, these could be
# house-keeping genes or a set of spike-in controls.

# If the negative control features are not known a priori, they can be identified
# empirically. This can be achieved via a 2-stage approach, RUVm.

# Stage 1 involves performing a differential methylation analysis using
# RUV-inverse (by default) and the 613 Illumina negative controls (INCs) as
# negative control features. This will produce a list of CpGs ranked by p-value
# according to their level of association with the factor of interest. This list
# can then be used to identify a set of empirical control probes (ECPs), which will
# capture more of the unwanted variation than using the INCs alone. ECPs are selected
# by designating a proportion of the CpGs least associated with the factor of interest
# as negative control features; this can be done based on either an FDR cut-off
# or by taking a fixed percentage of probes from the bottom of the ranked list.

#Stage 2 involves performing a second differential methylation analysis
# on the original data using RUV-inverse (by default) and the ECPs.

##### Starting Batch Normalisation #####

# setup the factor of interest
grp <- factor(targets$Sample_Group, labels=c(0,1))
# extract Illumina negative control data
INCs <- getINCs(rgSet)
# add negative control data to M-values
Mc <- rbind(Mval,INCs)
# create vector marking negative controls in data matrix
ctl1 <- rownames(Mc) %in% rownames(INCs)

# Stage 1 analysis - if no. of samples > 613 for 450k data
des <- model.matrix(~grp)
# limma differential methylation analysis
lfit1 <- lmFit(Mval, design=des)
lfit2 <- eBayes(lfit1) # Stage 1 analysis
top1 <- topTable(lfit2, num=Inf) # Removing intercept from test coefficients
ctl2 <- rownames(Mval) %in% rownames(top1[top1$adj.P.Val > 0.5,])


# Stage 2 analysis

# Perform RUV adjustment and fit
rfit3 <- RUVfit(Y = Mval, X = grp, ctl = ctl2)
rfit4 <- RUVadj(Y = Mval, fit = rfit3)

# To visualise the effect that the RUVm adjustment is having on the data,
# using an MDS plot for example, the getAdj function can be used to
# extract the adjusted values from the RUVm fit object produced by RUVfit.

Madj <- getAdj(Mval, rfit3) # get adjusted values

# The MDS plots below show how the relationship between the samples
# changes with and without RUVm adjustment. RUVm reduces the distance between
# the samples in each group by removing unwanted variation.

pdf(file = "/Users/sups/Downloads/R_Prog/COV/limma.pdf", width = 10, height = 10)
par(mfrow=c(1,2))
plotMDS(Mval, labels=targets$Sample_Name, col=as.integer(factor(targets$Sample_Group)),
        main="Unadjusted", gene.selection = "common")
legend("topleft",legend=c("Healthy","RA"),pch=16,cex=1,col=1:2)
plotMDS(Madj, labels=targets$Sample_Name, col=as.integer(factor(targets$Sample_Group)),
        main="Adjusted: RUV-inverse", gene.selection = "common")
legend("topleft",legend=c("Healthy","RA"),pch=16,cex=1,col=1:2)
dev.off()

######## TESTING FOR DIFFERENCES BETWEEN GROUP VARIATIONS. #########
# We may be interested in CpG sites that are consistently methylated in the
# normal samples, but variably methylated in the RA samples.

# CALLING FUNCTION "varFit": We are interested in testing for differential
# variability in the RA versus healthy group. Note that when we specify the
# coef parameter, which corresponds to the columns of the design matrix to be used
# for testing differential variability, we need to specify both the intercept and
# the second column. The ID variable is a nuisance parameter and not used when
# obtaining the absolute deviations, however it can be included in the linear
# modelling step. For methylation data, the function will take either a matrix of
# M-values, β values or a object as input. If β values are supplied, a logit transformation
# is performed. Note that as a default, varFit uses the robust setting in the
# limma framework, which requires the use of the statmod package.

fitvar <- varFit(Madj, design = design, coef = c(1,2))

table <- get_p_table(fitvar, coef = 2, number = dim(fitvar)[1])
write.csv(table, file = "/Users/sups/Downloads/R_Prog/COV/RawPval.csv")
