########################################################################
## This is an exploration example of genomicRanges 
########################################################################
#######   Data selected by te user #####################################
########################################################################
args <- commandArgs(trailingOnly = TRUE)
exp_mat_path <-args[1]
# exp_mat_path <-c("../Results/Matrices_splited_by_gene/ABCB1/TCGA_Basal_splited_by_the_expression_of_the_gene_ABCB1_25th_top_low.t
Data_path <- c("../Data_RNA-Analysis-tools")
setwd("/media/antonis/mount_me/Raul_Documents/github/Data_RNA-Analysis-tools")
################################################
######  Libraries needed    ####################
################################################
if (!require("BiocManager")) {
  install.packages("BiocManager")
  library(BiocManager)}

if (!require("GenomicRanges")) {
  BiocManager::install("GenomicRanges")
  library(GenomicRanges)}
if (!require("Rsamtools")) {
  BiocManager::install("Rsamtools")
  library(Rsamtools)}
if (!require("GenomicFeatures")) {
  BiocManager::install("GenomicFeatures")
  library(GenomicFeatures)}
if (!require("GenomicAlignments")) {
  BiocManager::install("GenomicAlignments")
  library(GenomicAlignments)}
if (!require("TxDb.Hsapiens.UCSC.hg19.knownGene")) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)}

################################################
######  Exercises           ####################
################################################

read <- GRanges( seqnames = c("19"),
                ranges = IRanges( start=c(44047464) ,
                                 end = c(44047499) ), strand =c("+"),
                seqlengths=c("19"=591289983)
                )
names(read) <-c("XRCC1")

read@ranges@NAMES
names(read)

read.df <- as.data.frame(read)
read.df
?GRanges

read.df$width
dim(read.df)


#getAnywhere(readGAlignments)
#h1b <- readGAlignments(paste0(Data_path,"/","hESC1_chr18.sam"))
h1b <- readGAlignments(paste0(Data_path,"/","hESC1_chr18.bam"))
str(h1b)
dim(h1b)
h1b[1:10,]
h1b[strand(h1b)=="+",]
class(h1b[strand(h1b)=="+",])