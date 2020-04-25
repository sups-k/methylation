# Bar graphs of hyper and hypo methylated genes

hyper <- read.csv(file = "/Users/sups/Downloads/R_Prog/RESULTS/DownReg_689/689_upmethyl_microarray.csv", header = TRUE)
hypo <- read.csv(file = "/Users/sups/Downloads/R_Prog/RESULTS/UpReg_689/689_downmethyl_microarray.csv", header = TRUE)

hyper_counts <- table(hyper$UCSC_RefGene_Name)
hypo_counts <- table(hypo$UCSC_RefGene_Name)

barplot(hypo_counts, las=2, cex.names = 0.8, ylim = c(0,5 + max(hypo_counts)))
title(ylab="Number of CpGs mapping to genes", line=2.5, cex.lab=0.9)
title(xlab="Genes", line=0, cex.lab=0.9)
title(main="Hypomethylated Genes", line = 0)

barplot(hyper_counts, las=2, cex.names = 0.8, ylim = c(0,5 + max(hyper_counts)))
title(ylab="Number of CpGs mapping to genes", line=2.5, cex.lab=0.9)
title(xlab="Genes", line=0, cex.lab=0.9)
title(main="Hypermethylated Genes", line = 0)
