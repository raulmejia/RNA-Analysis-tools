# This script receives a expression matrix and an associated annotation file and generates boxplots with dots over them
##########
# required libraries
##########
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("ggplot2")) {
  install.packages("ggplot2", ask =FALSE)
  library("ggplot2")
}
if (!require("reshape2")) {
  install.packages("reshape2", ask =FALSE)
  library("reshape2")
}
if (!require("utils")) {
  install.packages("utils", ask =FALSE)
  library("utils")
}
if (!require("limma")) {
  install.packages("limma", ask =FALSE)
  library("limma")
}


############
## Variables defined by the user
############
myargs <- commandArgs(trailingOnly = TRUE)

path_matrix <- myargs[1]
# path_matrix <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23_HNGC_collapsed_by_median.tsv"
path_matrix <- normalizePath( path_matrix )

path_annotations <-"/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23_HNGC_collapsed_by_median_annotation.tsv"

Label_for_yor_results <- myargs[2]    
# Label_for_yor_results <- "Run_19_01_2020_"

Folder_to_save_results <- myargs[3]
Folder_to_save_results <- normalizePath( Folder_to_save_results )
dir.create(Folder_to_save_results,recursive = TRUE)

annotation_file_path <- 
  
  
############
# Body's program
############
my_matrix <- read.table(path_matrix , header = TRUE )
my_matrix <- read.table(path_matrix , header = TRUE )

annot <- read.table(path_annotations , header = TRUE )

head(my_matrix)
plotDensities(log2(my_matrix+1), legend = FALSE) # antonis matrix have a nice shape 
mylog2mat <- log2(my_matrix+1) 


data.frame(as.vector(my_matrix["ABCA4",]),annot)
plotme <- my_matrix["ABCA4",]
colnames(plotme) <- NULL
rownames(plotme) <- NULL
cbind(my_matrix["ABCA4",], annot )

mydf<- as.data.frame(my_matrix)

log24plot <- annot 
as.numeric(mydf["ABCA4",])
annot$gene <- as.numeric(mydf["ABCA4",])
log24plot$gene <- as.numeric(mylog2mat["ABCA4",])
#colnames(annot)[3] <- "ABCA4"

ggplot( annot, aes(x=group , y = gene) )  +  geom_boxplot()
ggplot( log24plot, aes(x=group , y = gene) )  +  geom_boxplot()






annot <- data.frame(colnames(my_matrix),sub( '[0-9]*$', '', colnames(my_matrix) ))
colnames(annot) <- c("sample","group")
write.table(annot, file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23_HNGC_collapsed_by_median_annotation.tsv", sep="\t", row.names = FALSE)     



ggplot(plotdf, aes(x=Treatment, y=Lipid.Droplets, fill=Cell.line)) + 
  geom_boxplot() + ggtitle(Title) + theme_grey(base_size = 22)
ggplot(plotdf, aes(x=Cell.line, y=Lipid.Droplets, fill=Treatment)) + 
  geom_boxplot()+ ggtitle(Title) + theme_grey(base_size = 22)
dev.off()

plotdf <- mydf[, c(fixed_colnames,data_col_names[2])]
plotdf[,"Cell.line"] <- as.factor(plotdf[,"Cell.line"])
plotdf[,"Treatment"] <- as.factor(plotdf[,"Treatment"])
str(plotdf)
Title <- data_col_names[2]
pdf(file= paste0(Title,".pdf"))
ggplot(plotdf, aes(x=Treatment, y=Autophagosomes, fill=Cell.line)) + 
  geom_boxplot() + ggtitle(Title) + theme_grey(base_size = 22)
ggplot(plotdf, aes(x=Cell.line, y=Autophagosomes, fill=Treatment)) + 
  geom_boxplot()+ ggtitle(Title) + theme_grey(base_size = 22)
dev.off()

plotdf <- mydf[, c(fixed_colnames,data_col_names[3])]
plotdf[,"Cell.line"] <- as.factor(plotdf[,"Cell.line"])
plotdf[,"Treatment"] <- as.factor(plotdf[,"Treatment"])
str(plotdf)
Title <- data_col_names[3]
pdf(file= paste0(Title,".pdf"))
ggplot(plotdf, aes(x=Treatment, y=Lamellar.bodies, fill=Cell.line)) + 
  geom_boxplot() + ggtitle(Title) + theme_grey(base_size = 22)
ggplot(plotdf, aes(x=Cell.line, y=Lamellar.bodies, fill=Treatment)) + 
  geom_boxplot()+ ggtitle(Title) + theme_grey(base_size = 22)
dev.off()





