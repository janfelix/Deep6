#Read in Deep6 output file
out = read.delim("test_data.fasta_predict_250bp_deep6.txt")

#Assign group prediction for each sequence based the higehst score, the score also has to be 1.25 that of the median scroe for the sequence
library(matrixStats)
out$deep6 = ifelse(out$duplo>=rowMaxs(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")]))&out$duplo>1.25*rowMedians(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")])),"duplo",ifelse(out$euk>=rowMaxs(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")]))&out$euk>1.25*rowMedians(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")])),"euk",ifelse(out$mono>=rowMaxs(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")]))&out$mono>1.25*rowMedians(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")])),"mono",ifelse(out$pro>=rowMaxs(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")]))&out$pro>1.25*rowMedians(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")])),"pro", ifelse(out$vari>=rowMaxs(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")]))&out$vari>1.25*rowMedians(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")])),"vari", ifelse(out$ribo>=rowMaxs(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")]))&out$ribo>1.25*rowMedians(as.matrix(out[,c("duplo", "euk", "mono", "pro", "ribo", "vari")])),"ribo", "uncertain"))))))
out$deep6 = factor(out$deep6, levels=c("euk", "pro", "duplo", "vari", "mono", "ribo", "uncertain"))

#Overview of contig predictions
table(out$deep6)