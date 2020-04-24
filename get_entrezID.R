# Loading all required libraries
library(missMethyl)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(org.Hs.eg.db)

# Read the table containing CpG probes with p-values
table <- read.csv("/Users/sups/Downloads/R_Prog/RawPval_GSE111942.csv", header = TRUE)
rownames(table) <- table$Probe_ID

# The following code is adapted from .getFlatAnnotation() function of the missMethyl package
# Author of .getFlatAnnotation() is Jovana Maksimovic
# Author of code below: Suparna Kumar

sig.cpg = rownames(table)
array.type = "450K"
# array.type = "EPIC"
sig.cpg <- as.character(sig.cpg)
sig.cpg <- sig.cpg[!is.na(sig.cpg)]

anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann.keep<-anno[grepl("^cg",anno$Name),]
missing<-ann.keep$UCSC_RefGene_Name==""
ann.keep<-ann.keep[!missing,]

geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
names(geneslist)<-rownames(ann.keep)
grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
names(grouplist)<-rownames(ann.keep)
# We are not getting every CpG mapped to the genes in UCSC Browser
flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))

flat$symbol<-as.character(flat$symbol)
flat$group <- as.character(flat$group)
flat$cpg<- substr(rownames(flat),1,10)
flat$alias <- suppressWarnings(limma:::alias2SymbolTable(flat$symbol))

eg <- suppressMessages(select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c("ENTREZID","SYMBOL"), keytype="ENTREZID"))
colnames(eg) <- c("gene_id","symbol")

flat$cpg<- rownames(flat)

m <- match(flat$alias,eg$symbol)
flat$entrezid <- eg$gene_id[m]

flat <- flat[!is.na(flat$entrezid),]
id<-paste(flat$cpg,flat$entrezid,sep=".")
d <- duplicated(id)
flat.u <- flat[!d,]

# Annotate Entrez IDs to the table of p-values
a <- merge(table, flat.u, by = 0) # merge by row names
a <- a[, -11]
rownames(a) <- a$Row.names
a <- a[, -1]
a <- a[, -1] # This is not a mistake. I want to remove columns 1 & 2. After removing column 1, the column 2 becomes the new column 1
names(a)[names(a) == "symbol"] <- "gene_name"

# Write to CSV file
write.csv(a, file = "/Users/sups/Downloads/R_Prog/RawPval_GSE111942_AnnotatedEntrez.csv")

######## From manifest file #########
anno <- read.csv("/Users/sups/Downloads/R_Prog/450k_gene_manifest.csv", header = TRUE)
# then read your p value file as the variable "table"
anno.keep <- anno[grepl("^cg", anno$Name),]
anno.keep <- anno.keep[,-1]
rownames(anno.keep) <- anno.keep$Name

a <- merge(table, anno.keep, by = 0)

names(a)[names(a) == "Row.names"] <- "CpG"

# To remove any column, use a <- a[, -5] where the col no. = 5
sigCpG <- a[a$p.BH_X1.1 < 0.05,] # all significant CpG
upCpG <- sigCpG[sigCpG$mean > 0,] # hypermethylated
downCpG <- sigCpG[sigCpG$mean < 0,] # hypomethylated
upCpG <- upCpG[order(upCpG[,10], decreasing = TRUE),] # sort in descending order of M-values
downCpG <- downCpG[order(downCpG[,10], decreasing = TRUE),]

write.csv(upCpG, file = "/Users/sups/Downloads/R_Prog/689_filtered_hyPERmethylation.csv")
write.csv(downCpG, file = "/Users/sups/Downloads/R_Prog/689_filtered_hyPOmethylation.csv")
write.csv(a, file = "/Users/sups/Downloads/R_Prog/Pval_GSE42861_AnnotatedUCSC.csv")