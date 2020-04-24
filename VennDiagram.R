# Venn diagram to depict the common genes and probes between the datasets
# GSE42861 and GSE111942

# Load library
library(VennDiagram)

# Load the annotated CpG result tables - change to hypermethylation if required

# res_689 <- read.csv("/Users/sups/Downloads/R_Prog/Pval_GSE42861_AnnotatedUCSC.csv")
# res_43 <- read.csv("/Users/sups/Downloads/R_Prog/Pval_GSE111942_AnnotatedUCSC.csv")

res_689 <- read.csv("/Users/sups/Downloads/R_Prog/RESULTS/results_689/689_filtered_hyPOmethylation.csv")
res_43 <- read.csv("/Users/sups/Downloads/R_Prog/RESULTS/results_43/43_filtered_hyPOmethylation.csv")
common_CpG <- read.csv("/Users/sups/Downloads/R_Prog/RESULTS/689_43_common_CpG_hyPOmethylated.csv")

# Remove the column named "X" which contains row indices - check if you have it otherwise
# DO NOT run the code below
res_689 <- res_689[,-1]
res_43 <- res_43[,-1]
common_CpG <- common_CpG[,-1]

# # Get genes with p.BH < 0.05 and mean < -1 ==> hypomethylated
# res_689 <- res_689[res_689$p.BH_X1.1 < 0.05,]
# res_43 <- res_43[res_43$p.BH_X1.1 < 0.05,]
# 
# res_689 <- res_689[res_689$mean < 0,]
# res_43 <- res_43[res_43$mean < 0,]

# Sort gene names in ascending order
res_689 <- res_689[order(res_689[,"UCSC_RefGene_Name"]),]
res_43 <- res_43[order(res_43[,"UCSC_RefGene_Name"]),]
common_CpG <- common_CpG[order(common_CpG[,"UCSC_RefGene_Name"]),]

# Convert gene names to character and select them
gene_689 = as.character(unique(res_689$UCSC_RefGene_Name))
gene_43 = as.character(unique(res_43$UCSC_RefGene_Name))
common_gene = as.character(unique(common_CpG$UCSC_RefGene_Name))

# Remove blank space
gene_689 <- gene_689[-1]
gene_43 <- gene_43[-1]
common_gene <- common_gene[-1]

# Write as text files
write(gene_689, file = "/Users/sups/Downloads/R_Prog/gene_689.txt")
write(gene_43, file = "/Users/sups/Downloads/R_Prog/gene_43.txt")
write(common_gene, file = "/Users/sups/Downloads/R_Prog/683_49_hypo.txt")

# Convert CpG probe names to character and select them
probe_689 = as.character(unique(res_689$CpG))
probe_43 = as.character(unique(res_43$CpG))
probe_gene = as.character(unique(common_CpG$CpG))

# Remove some variables
rm(common_gene, gene_689, gene_43)

# Linux commands to remove semicolons
system(paste("tr ';' '\n' < /Users/sups/Downloads/R_Prog/gene_689.txt > /Users/sups/Downloads/R_Prog/out_gene_689.txt"))
system(paste("tr ';' '\n' < /Users/sups/Downloads/R_Prog/gene_43.txt > /Users/sups/Downloads/R_Prog/out_gene_43.txt"))
system(paste("tr ';' '\n' < /Users/sups/Downloads/R_Prog/683_49_hypo.txt > /Users/sups/Downloads/R_Prog/out_683_49_hypo.txt"))

# Remove old files and rename the new ones
system(paste("rm /Users/sups/Downloads/R_Prog/gene_689.txt"))
system(paste("rm /Users/sups/Downloads/R_Prog/gene_43.txt"))
system(paste("rm /Users/sups/Downloads/R_Prog/683_49_hypo.txt"))

system(paste("mv /Users/sups/Downloads/R_Prog/out_gene_689.txt /Users/sups/Downloads/R_Prog/gene_689.txt"))
system(paste("mv /Users/sups/Downloads/R_Prog/out_gene_43.txt /Users/sups/Downloads/R_Prog/gene_43.txt"))
system(paste("mv /Users/sups/Downloads/R_Prog/out_683_49_hypo.txt /Users/sups/Downloads/R_Prog/683_49_hypo.txt"))

# Sort and filter out unique gene names
system(paste("sort /Users/sups/Downloads/R_Prog/gene_689.txt | uniq > /Users/sups/Downloads/R_Prog/unique_gene_689.txt"))
system(paste("sort /Users/sups/Downloads/R_Prog/gene_43.txt | uniq > /Users/sups/Downloads/R_Prog/unique_gene_43.txt"))
system(paste("sort /Users/sups/Downloads/R_Prog/683_49_hypo.txt | uniq > /Users/sups/Downloads/R_Prog/unique_683_49_hypo.txt"))

# Load the gene names
names_689 <- read.delim("/Users/sups/Downloads/R_Prog/unique_gene_689.txt", header = FALSE)
names_43 <- read.delim("/Users/sups/Downloads/R_Prog/unique_gene_43.txt", header = FALSE)
names_hypo <- read.delim("/Users/sups/Downloads/R_Prog/unique_683_49_hypo.txt", header = FALSE)

# Convert gene names to character
colnames(names_43) <- "Gene"
colnames(names_689) <- "Gene"
colnames(names_hypo) <- "Gene"
names_689 <- as.character(names_689$Gene)
names_43 <- as.character(names_43$Gene)
names_hypo <- as.character(names_hypo$Gene)

# Venn diagrams of gene names

# Scaled Venn diagram with numbers
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Scaled_Hypomethylated_Genes_num.png", width = 500, height = 500)
draw.pairwise.venn(length(names_689), length(names_43), length(names_hypo), category = c("GSE42861 Genes", "GSE111942 Genes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
dev.off()

# Scaled Venn diagram with percentage
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Scaled_Hypomethylated_Genes_percent.png", width = 500, height = 500)
draw.pairwise.venn(length(names_689), length(names_43), length(names_hypo), category = c("GSE42861 Genes", "GSE111942 Genes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = TRUE, print.mode = "percent")
dev.off()

# Unscaled Venn diagram with percentage
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Hypomethylated_Genes_num.png", width = 500, height = 500)
draw.pairwise.venn(length(names_689), length(names_43), length(names_hypo), category = c("GSE42861 Genes", "GSE111942 Genes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = FALSE, print.mode = "percent")
dev.off()

# Unscaled Venn diagram with numbers
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Hypomethylated_Genes_percent.png", width = 500, height = 500)
draw.pairwise.venn(length(names_689), length(names_43), length(names_hypo), category = c("GSE42861 Genes", "GSE111942 Genes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = FALSE)
dev.off()

# Venn diagrams of CpG probes

# Scaled Venn diagram with numbers
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Scaled_Hypomethylated_CpG_num.png", width = 500, height = 500)
draw.pairwise.venn(length(probe_689), length(probe_43), length(probe_gene), category = c("GSE42861 CpG Probes", "GSE111942 CpG Probes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = TRUE)
dev.off()

# Scaled Venn diagram with percentage
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Scaled_Hypomethylated_CpG_percent.png", width = 500, height = 500)
draw.pairwise.venn(length(probe_689), length(probe_43), length(probe_gene), category = c("GSE42861 CpG Probes", "GSE111942 CpG Probes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = TRUE, print.mode = "percent")
dev.off()

# Unscaled Venn diagram with numbers
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Hypomethylated_CpG_num.png", width = 500, height = 500)
draw.pairwise.venn(length(probe_689), length(probe_43), length(probe_gene), category = c("GSE42861 CpG Probes", "GSE111942 CpG Probes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = FALSE)
dev.off()

# Unscaled Venn diagram with percentage
png(filename = "/Users/sups/Downloads/R_Prog/RESULTS/Hypomethylated_CpG_percent.png", width = 500, height = 500)
draw.pairwise.venn(length(probe_689), length(probe_43), length(probe_gene), category = c("GSE42861 CpG Probes", "GSE111942 CpG Probes"), lty = rep("blank", 2), fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2), scaled = FALSE, print.mode = "percent")
dev.off()