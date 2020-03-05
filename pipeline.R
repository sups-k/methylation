# References:
# https://hhhh5.github.io/ewastools/articles/exemplary_ewas.html
# https://rdrr.io/github/hhhh5/ewastools/man/detectionP.html
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6346546/

rm(list=ls())

# Loading all required libraries
library(missMethyl)
library(limma)
library(minfi)
library(ewastools)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(data.table)

# Set directory containing all IDAT files and Sample Sheet
baseDir <- "/Users/sups/Documents/R_Prog/COV/GSE42861"
setwd(baseDir)

# Read the Sample Sheet
# DOES NOT REQUIRE IDAT FILES IN SENTRIX_ID FOLDERS
targets <- read.metharray.sheet(baseDir)

# Display all columns of Sample Sheet except the Sentrix ID, Sentrix Position, file path
targets[,1:5]

# Display Sentrix ID, Sentrix Position, file path
targets[,6:8]

##################################################################################
## Calculate detection P values - adapted from ewastools
# Author: Jonathan A. Heiss

# Reading the files into object "meth" using ewastools
# For this, all the files must be in the same folder, not separated as in minfi
meth = read_idats(targets$Basename,quiet=TRUE)

# Calculate detection p-values for each locus
meth = detectionP(meth)

# Define number of males and females in the sample
males =  which(targets$sex=="m")
females =  which(targets$sex=="f")

# Warnings defined in case sex or detP information is missing
if(is.null(males) | is.null(females)) stop('Please specify the column indices for male and female subjects')
if(!'detP'%in%names(meth)) stop('detP component missing')

# Creating table of Y chromosome probe intensities from the data
chrY1 = meth$manifest[chr=='Y',index]
probe <- meth$manifest$probe_id
probe_position=which(probe=="cg13618458")
chrY=chrY1[!chrY1 %in% probe_position]
chrY = meth$detP[chrY,]

# Some general p-value cut-offs used in methylation
cutoffs = c(1,0.5,0.1,0.05,0.01,0.001,0.0001)

# Calculating quantile values for males and females to find out how many Y chromosome
# probes are undetected in both sexes after applying the respective p-value cut-off.
# For females, this number should be as high as possible and for males it should be very low
tmp = sapply(cutoffs,function(t){ colSums(chrY>t,na.rm=TRUE) })
males = apply(tmp[males  ,],2,quantile,prob=0.9)
females = apply(tmp[females,],2,quantile,prob=0.1)

# Plotting the above results
plot(-log10(cutoffs),females,ylim=c(0,nrow(chrY)),ylab='Chr Y # undetected ',xlab='p-value cutoff',xaxt="n")
points(-log10(cutoffs),males,pch=3)
axis(1,at=-log10(cutoffs),labels=cutoffs)
legend('topleft',pch=c(2,1),legend=c('Male 90% Quantile','Female 10% Quantile'))
invisible(NULL)

#############################################################################################
# CHOOSE YOUR P-VALUE #

# Remove the large objects that won't be used
rm(list = c("meth", "males", "females", "chrY", "cutoffs", "tmp"))
detach("package:ewastools", unload = TRUE)

# Import the IDAT data into an RGChannel object for further analysis
rgSet <- read.metharray.exp(targets = targets)

# Converting to MethylSet object & performing SWAN normalization
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=FALSE)

# Plotting density distribution of beta values before and after using SWAN.
par(mfrow=c(1,1), cex=0.50)
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")

# Only select those intensities that are lesser than the cut-off
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.05) == ncol(rgSet) # where 0.001 is the selected p-value
mSetSw <- mSetSw[keep,]

# Extract beta and M-values from the SWAN normalised data.
# We prefer to add an offset to the methylated and unmethylated intensities
# when calculating M-values, hence we extract the methylated and unmethylated channels
# separately and perform our own calculation.
meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSw)
dim(Mval)

# Plot MDS (multi-dimensional scaling) of cancer and normal samples.
# This is a good check to make sure samples cluster together according to their type.
par(mfrow=c(1,1))
plotMDS(Mval, labels=targets$Sample_Name, col=as.integer(factor(targets$Sample_Group)))
legend("bottom",legend=c("Healthy","RA"),pch=16,cex=1.2,col=1:2)

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
par(mfrow=c(1,1))
for(i in 1:4){
  stripchart(beta[rownames(beta)==cpgs[i],]~design[,2],method="jitter",
             group.names=c("Healthy","RA"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
             vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
  title(cpgs[i],cex.main=1.5)
}

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

# setup the factor of interest
grp <- factor(targets$Sample_Group, labels=c(0,1))
# extract Illumina negative control data
INCs <- getINCs(rgSet)
head(INCs)
# add negative control data to M-values
Mc <- rbind(Mval,INCs)
# create vector marking negative controls in data matrix
ctl1 <- rownames(Mc) %in% rownames(INCs)
table(ctl1)

# Stage 1 analysis
rfit1 <- RUVfit(Y = Mc, X = grp, ctl = ctl1)
rfit2 <- RUVadj(Y = Mc, fit = rfit1)

top1 <- topRUV(rfit2, num=Inf, p.BH = 1) # p.BH is cutoff value for Benjamini-Hochberg adjusted p-values
head(top1)

ctl2 <- rownames(Mval) %in% rownames(top1[top1$p.BH_X1.1 > 0.5,])
table(ctl2)

# Stage 2 analysis

# Perform RUV adjustment and fit
rfit3 <- RUVfit(Y = Mval, X = grp, ctl = ctl2)
rfit4 <- RUVadj(Y = Mval, fit = rfit3)
# Look at table of top results
topRUV(rfit4)

# To visualise the effect that the RUVm adjustment is having on the data,
# using an MDS plot for example, the getAdj function can be used to
# extract the adjusted values from the RUVm fit object produced by RUVfit.
# NOTE: The adjusted values should only be used for visualisations - it is
# NOT recommended that they are used in any downstream analysis.

Madj <- getAdj(Mval, rfit3) # get adjusted values

# The MDS plots below show how the relationship between the samples
# changes with and without RUVm adjustment. RUVm reduces the distance between
# the samples in each group by removing unwanted variation.

par(mfrow=c(1,2))
plotMDS(Mval, labels=targets$Sample_Name, col=as.integer(factor(targets$Sample_Group)),
        main="Unadjusted", gene.selection = "common")
legend("topleft",legend=c("RA","Healthy"),pch=16,cex=1,col=1:2)
plotMDS(Madj, labels=targets$Sample_Name, col=as.integer(factor(targets$Sample_Group)),
        main="Adjusted: RUV-inverse", gene.selection = "common")
legend("topleft",legend=c("RA","Healthy"),pch=16,cex=1,col=1:2)

# Rather than testing for differences in mean methylation, we may be interested
# in testing for differences between group variances. For example, it has been
#hypothesised that highly variable CpGs in cancer are important for tumour
# progression (Hansen et al. 2011). Hence we may be interested in CpG sites that
# are consistently methylated in the normal samples, but variably methylated in the
# cancer samples.

# In general we recommend at least 10 samples in each group for accurate variance
# estimation. We are interested in testing for differential
# variability in the RA versus healthy group. Note that when we specify the
# coef parameter, which corresponds to the columns of the design matrix to be used
# for testing differential variability, we need to specify both the intercept and
# the second column. The ID variable is a nuisance parameter and not used when
# obtaining the absolute deviations, however it can be included in the linear
# modelling step. For methylation data, the function will take either a matrix of
# M-values, β values or a object as input. If β values are supplied, a logit transformation
# is performed. Note that as a default, varFit uses the robust setting in the
# limma framework, which requires the use of the statmod package.

fitvar <- varFit(Mval, design = design, coef = c(1,2))

# The numbers of hyper-variable (1) and hypo-variable (-1) genes in RA vs healthy
# can be obtained using decideTests.

summary(decideTests(fitvar))
topDV <- topVar(fitvar, coef=2)
topDV

# The β values for the top 4 differentially variable CpGs can be seen below:
cpgsDV <- rownames(topDV)
par(mfrow=c(2,2))
for(i in 1:4){
  stripchart(beta[rownames(beta)==cpgsDV[i],]~design[,2],method="jitter",
             group.names=c("Healthy","RA"),pch=16,cex=1.5,col=c(4,2),ylab="Beta values",
             vertical=TRUE,cex.axis=1.5,cex.lab=1.5)
  title(cpgsDV[i],cex.main=1.5)
}

## Gene Ontology Analysis

# Once a differential methylation analysis has been performed, it may be
# of interest to know which gene pathways are targeted by the significant
# CpG sites.

# For the Illumina Infinium HumanMethylation450 array the number of probes
# per gene ranges from 1 to 1299, with a median of 15 probes per gene. For 
# the EPIC array, the range is 1 to 1487, with a median of 20 probes per gene.
# This means that when mapping CpG sites to genes, a gene is more likely to be
# selected if there are many CpG sites associated with the gene.

# One way to take into account this selection bias is to model the relationship
# between the number of probes per gene and the probability of being selected.
# This can be performed by adapting the goseq method of Young et al. Each gene
# then has a prior probability associated with it, and a modified version of a
# hypergeometric test can be performed, testing for over-representation of the
# selected genes in each gene set.

# The gometh function performs gene set testing on GO categories or KEGG pathways
# The gsameth function is a more generalised gene set testing function which can
# take as input a list of user specified gene sets. Note that for gsameth, the
# format for the gene ids for each gene in the gene set needs to be Entrez Gene
# IDs. For example, the entire curated gene set list (C2) from the Broad’s
# Molecular Signatures Database can be specified as input. The R version of these
# lists can be downloaded from http://bioinf.wehi.edu.au/software/MSigDB/index.html.
# Both functions take a vector of significant CpG probe names as input.

# To illustrate how to use gometh, consider the results from the differential
# methylation analysis with RUVm.
top <- topRUV(rfit4, number = Inf, p.BH = 1)
table(top$p.BH_X1.1 < 0.01)

# We take the top 10000 CpG sites as input to `gometh`.
topCpGs<-topRUV(rfit4,number=10000)
sigCpGs <- rownames(topCpGs)

# The takes as input a character vector of CpG names, and optionally, a
# character vector of all CpG sites tested. If the all.cpg argument is
# omitted, all the CpGs on the array are used as background. To change the
# array type, the array.type argument can be specified as either “450K” or
# “EPIC”. The default is “450K”.

# If the plot.bias argument is TRUE, a figure showing the relationship between
# the probability of being selected and the number of probes per gene will be
# displayed.

# For testing of GO terms, the collection argument takes the value “GO”, which
# is the default setting. For KEGG pathway analysis, set collection to “KEGG”.
# The function topGSA shows the top enriched GO categories. The function gsameth
# is called for GO and KEGG pathway analysis with the appropriate inputs.

gst <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(top), collection="GO")
gene_ontology_res <- topGSA(gst) # top 20 pathways
write.csv(gene_ontology_res, file = "Gene_Ontology_GSE42861-60.csv")
entrez_ids <- getMappedEntrezIDs(sig.cpg = sigCpGs, all.cpg = rownames(top), array.type = "450K")
# It maps the significant CpG probe names to Entrez Gene IDs, as well as all the CpG sites tested.
# It also calculates the numbers of probes for gene.