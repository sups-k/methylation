table1 <- read.csv("/Users/sups/Downloads/R_Prog/GSE45291_hgu133plus2_RA_control_annotated.csv", header = TRUE)
table2 <- read.csv("/Users/sups/Downloads/R_Prog/AdjPval_GSE48261_AnnotatedEntrez.csv", header = TRUE)
table3 <- read.csv("/Users/sups/Downloads/R_Prog/AdjPval_GSE111942_AnnotatedEntrez.csv", header = TRUE)

a <- table1
b <- table2
c <- table3

a <- a[which(a$FDR_corrected < 0.05),]
b <- b[which(b$Adj.P.Value < 0.05),]
c <- c[which(c$Adj.P.Value < 0.05),]

b1 <- b[which(b$LogVarRatio > 1),]
b2 <- b[which(b$LogVarRatio < -1),]

c1 <- c[which(c$LogVarRatio > 1),]
c2 <- c[which(c$LogVarRatio < -1),]

a1 <- a[which(a$FC_log2 > 1),]
a2 <- a[which(a$FC_log2 < -1),]

a1 <- a1[order(a1[,5], decreasing = TRUE),]
b1 <- b1[order(b1[,3], decreasing = TRUE),]
c1 <- c1[order(c1[,3], decreasing = TRUE),]

a2 <- a2[order(a2[,5]),]
b2 <- b2[order(b2[,3]),]
c2 <- c2[order(c2[,3]),]

write.csv(a1, file = "/Users/sups/Downloads/R_Prog/revathi_filtered_up.csv")
write.csv(a2, file = "/Users/sups/Downloads/R_Prog/revathi_filtered_down.csv")

write.csv(b1, file = "/Users/sups/Downloads/R_Prog/689_filtered_up.csv")
write.csv(b2, file = "/Users/sups/Downloads/R_Prog/689_filtered_down.csv")

write.csv(c1, file = "/Users/sups/Downloads/R_Prog/43_filtered_up.csv")
write.csv(c2, file = "/Users/sups/Downloads/R_Prog/43_filtered_down.csv")