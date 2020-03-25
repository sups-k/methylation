## Calculate detection P values - adapted from ewastools
# Author: Jonathan A. Heiss

library(minfi)
library(data.table)
library(ewastools)

# Set directory containing all IDAT files and Sample Sheet
baseDir <- "/Users/sups/Downloads/R_Prog/COV/GSE42861"
setwd(baseDir)

# Read the Sample Sheet
targets <- read.metharray.sheet(baseDir)

# Reading the files into object "meth" using ewastools
meth = read_idats(targets$Basename,quiet=TRUE)

# Calculate detection p-values for each locus
meth = ewastools::detectionP(meth)

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
pdf(file = "/Users/sups/Documents/R_Prog/COV/DetP.pdf", width = 10, height = 10)

plot(-log10(cutoffs),females,ylim=c(0,nrow(chrY)),ylab='Chr Y # undetected ',xlab='p-value cutoff',xaxt="n")
points(-log10(cutoffs),males,pch=3)
axis(1,at=-log10(cutoffs),labels=cutoffs)
legend('topleft',pch=c(2,1),legend=c('Male 90% Quantile','Female 10% Quantile'))
invisible(NULL)

dev.off()