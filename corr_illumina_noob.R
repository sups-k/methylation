rm(list=ls())

# Loading all required libraries
library(missMethyl)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Set directory containing all IDAT files and Sample Sheet
baseDir1 <- "/Users/sups/Downloads/R_Prog/COV/GSE42861"
baseDir2 <- "/Users/sups/Downloads/R_Prog/COV/GSE"

# Read the Sample Sheet
targets1 <- read.metharray.sheet(baseDir1)
targets2 <- read.metharray.sheet(baseDir2)

#### Step 1: Read the raw IDAT files #####

# Import the IDAT data into an RGChannel object
rgSet1 <- read.metharray.exp(targets = targets1)
rgSet2 <- read.metharray.exp(targets = targets2)

#### Step 2: Convert the raw intensities into beta values #####

# Converting to MethylSet object
mSet1 <- preprocessIllumina(rgSet1) # Illumina background correction for rgSet1
mSet2 <- preprocessIllumina(rgSet2) # Illumina background correction for rgSet2

mSet3 <- preprocessNoob(rgSet1) # Noob background correction for rgSet1
mSet4 <- preprocessNoob(rgSet2) # Noob background correction for rgSet1

# Extract beta and M-values from the MethylSet data.
# We prefer to add an offset to the methylated and unmethylated intensities
# when calculating M-values, hence we extract the methylated and unmethylated channels
# separately and perform our own calculation.

# rgSet1 Illumina
meth1 <- getMeth(mSet1)
unmeth1 <- getUnmeth(mSet1)
colnames(meth1) <- sub("\\_.*", "", colnames(meth1))
colnames(unmeth1) <- sub("\\_.*", "", colnames(unmeth1))

# rgSet2 Illumina
meth2 <- getMeth(mSet2)
unmeth2 <- getUnmeth(mSet2)
colnames(meth2) <- sub("\\_.*", "", colnames(meth2))
colnames(unmeth2) <- sub("\\_.*", "", colnames(unmeth2))

# rgSet1 Noob
meth3 <- getMeth(mSet3)
unmeth3 <- getUnmeth(mSet3)
colnames(meth3) <- sub("\\_.*", "", colnames(meth3))
colnames(unmeth3) <- sub("\\_.*", "", colnames(unmeth3))

# rgSet2 Noob
meth4 <- getMeth(mSet4)
unmeth4 <- getUnmeth(mSet4)
colnames(meth4) <- sub("\\_.*", "", colnames(meth4))
colnames(unmeth4) <- sub("\\_.*", "", colnames(unmeth4))

Mval1 <- log2((meth1 + 100)/(unmeth1 + 100)) # rgSet1 Illumina
Mval2 <- log2((meth2 + 100)/(unmeth2 + 100)) # rgSet2 Illumina
Mval3 <- log2((meth3 + 100)/(unmeth3 + 100)) # rgSet1 Noob
Mval4 <- log2((meth4 + 100)/(unmeth4 + 100)) # rgSet2 Noob

#### Plotting the correlation ####

# rgSet1 - Plot correlation and save to PDF
pdf(file = "rgSet1_correlation.pdf", width = 10, height = 10)
for (i in 1:ncol(Mval1)) {
  plot(Mval1[, i], Mval3[, i], main="Scatterplot Noob vs Illumina", xlab="Illumina Mval", ylab="Noob Mval", pch=19)
  abline(lm(Mval3[, i]~Mval1[, i]), col="red")
}

# rgSet2 - Plot correlation and save to PDF
pdf(file = "rgSet2_correlation.pdf", width = 10, height = 10)
for (i in 1:ncol(Mval2)) {
  plot(Mval2[, i], Mval4[, i], main="Scatterplot Noob vs Illumina", xlab="Illumina Mval", ylab="Noob Mval", pch=19)
  abline(lm(Mval4[, i]~Mval2[, i]), col="red")
}