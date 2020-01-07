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
if (!require("SNPlocs.Hsapiens.dbSNP.20120608")) {
  BiocManager::install("SNPlocs.Hsapiens.dbSNP.20120608")
  library(SNPlocs.Hsapiens.dbSNP.20120608)}
if (!require("BSgenome")) {
  BiocManager::install("BSgenome")
  library(BSgenome)}
if (!require("BSgenome.Hsapiens.UCSC.hg19")) {
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  library(BSgenome.Hsapiens.UCSC.hg19)}

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


txdb<- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene ,
                     "gene")
str(txdb)
class(txdb)
txdb
exonsdb<- transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene ,
                        "exon")


hits <-  countOverlaps(h1b, txdb)
ol <- countOverlaps( txdb , h1b[hits==1])

class(ol)
ol
table(ol)

# Representing SNPs on R
                                        # The SNPlocs.Hsapiens.dbSNP.20120608 package is deprecated. Please use a SNPlocs data package based on a more recent dbSNP BUILD instead (e.g. BUILD 144 or BUILD 149). You can call BSgenome::available.SNPs()

genome <-  injectSNPs( BSgenome.Hsapiens.UCSC.hg19 , "SNPlocs.Hsapiens.dbSNP.20120608")
seqnames(genome)

