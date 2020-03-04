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
baseDir <- "/Users/sups/Documents/R_Prog/COV/Meth_Data"
setwd(baseDir)

# Read the Sample Sheet
# DOES NOT REQUIRE IDAT FILES IN SENTRIX_ID FOLDERS
targets <- read.metharray.sheet(baseDir)

# Display all columns of Sample Sheet except the Sentrix ID, Sentrix Position, file path
targets[,1:5]

# Display Sentrix ID, Sentrix Position, file path
targets[,6:8]

# Import the IDAT data into an RGChannel object for further analysis
rgSet <- read.metharray.exp(targets = targets)

# Converting to MethylSet object & performing SWAN normalization
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=FALSE)

# Plotting density distribution of beta values before and after using SWAN.
par(mfrow=c(1,1), cex=0.50)
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")

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
legend('topleft',pch=c(3,1),legend=c('Male 90% Quantile','Female 10% Quantile'))
invisible(NULL)

#############################################################################################
# CHOOSE YOUR P-VALUE #

# Remove the large objects that won't be used
rm(list = c("meth", "males", "females", "chrY", "cutoffs", "tmp"))

# Only select those intensities that are lesser than the cut-off
detP <- detectionP.minfi(rgSet)
keep <- rowSums(detP < 0.001) == ncol(rgSet) # where 0.001 is the selected p-value
mSetSw <- mSetSw[keep,]

# Extract beta and M-values from the SWAN normalised data.
# We prefer to add an offset to the methylated and unmethylated intensities
# when calculating M-values, hence we extract the methylated and unmethylated channels
# separately and perform our own calculation. For all subsequent analysis we use a random
# selection of 20000 CpGs to reduce computation time.
mset_reduced <- mSetSw[sample(1:nrow(mSetSw), 20000),]
meth <- getMeth(mset_reduced)
# meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mset_reduced)
# unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mset_reduced)
# beta <- getBeta(mSetSw)
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
## hello
