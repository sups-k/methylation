# Quality control - sex prediction
# Developed by Suparna Kumar
# This is to make sure that the predicted sex matches the reported sex for the data
# In case of a MISMATCH, please USE the PREDICTED SEX value

# Martin J. Aryee, Andrew E. Jaffe, Hector Corrada-Bravo, Christine Ladd-Acosta,
# Andrew P. Feinberg, Kasper D. Hansen, and Rafael A. Irizarry, Minfi: A flexible
# and comprehensive Bio- conductor package for the analysis of Infinium DNA Methylation
# microarrays, Bioinformatics 30 (2014), no. 10, 1363–1369.

rm(list=ls())

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Set directory containing all IDAT files and Sample Sheet
baseDir <- "/Users/sups/Downloads/R_Prog/COV/GSE42861"
setwd(baseDir)

# Read the Sample Sheet
targets <- read.metharray.sheet(baseDir)

# Import the IDAT data into an RGChannel object
rgSet <- read.metharray.exp(targets = targets)

mSex <- preprocessRaw(rgSet)
GRsex <- mapToGenome(mSex)
predictedSex <- getSex(GRsex, cutoff = -2)$predictedSex
sex <- cbind.data.frame(targets$sex, predictedSex, row.names = targets$Sample_Name)
write.csv(sex, file = "/Users/sups/Downloads/R_Prog/COV/predictedSex.csv")