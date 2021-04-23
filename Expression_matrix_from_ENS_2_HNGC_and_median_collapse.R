# This script converts the ENSBL to genesymbols from a matrix 
# The final numbers after the point are the version
# https://m.ensembl.org/Help/Faq?id=488

## Notes
## this program erase the numbers after the point in ENS codes this left the .XXX_PAR_region
## untouched if you want to erase it you can substitute for this command ensemblids_noversion <- sub( '\\.[0-9]*$', '', ensemblids_version )
## for this one  ensemblids_noversion <- sub( '\\.[[:graph:]]*$', '', ensemblids_version ) # Deletes PAR label but creates replicated 
## Nevertheless an undesirable effect of this chanve is the duplicity in ENS codes before the point

## Pending tasks, in the "filling matrix section" you can move it to a separated function and/or improve it with parallel 
## functions because is the most time consuming section
############################## 
## Required libraries
##############################
if (!require("gprofiler2")) {
  BiocManager::install("gprofiler2", dependencies = TRUE)
  library("gprofiler2")
}
if (!require("argparse")) {
  BiocManager::install("argparse", dependencies = TRUE)
  library("argparse")
}
if (!require("robustbase")) {
  BiocManager::install("robustbase", dependencies = TRUE)
  library("robustbase")
}
if (!require("dplyr")) {
  BiocManager::install("dplyr", dependencies = TRUE)
  library("dplyr")
}



############################## 
## Data given by the user
##############################
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-q", "--quietly", action="store_false", 
                    dest="verbose", help="Print little output")
parser$add_argument("-i", "--inputfile", type="character", 
                    help="input file with your gene list in genesymbols")
parser$add_argument("-o", "--outputfile", type="character", 
                    help="output file where you want to store your results")
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()
# print some progress messages to stderr if "quietly" wasn't requested

#############################
## The program starts
#############################
inputmatrix <-read.table( file=args$inputfile, stringsAsFactors = FALSE )
# inputmatrix <-read.table(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23.txt", stringsAsFactors = FALSE)

path2save <- args$outputfile
#  path2save <- "/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23_HNGC_collapsed_by_median.tsv"

ensemblids_version <- rownames( inputmatrix )
# ensemblids_noversion <- sub( '\\.[0-9]*$', '', ensemblids_version )

# getting the map between ESN (ensembl) and HGNC
table_HGNCids <- gprofiler2::gconvert( ensemblids_noversion , organism = "hsapiens", target="HGNC" )

# Checking the table with the converted Ids
# Identifying the input ids that mapped to more than 1 HNGC
duplicatedinputids <- table_HGNCids$input[ duplicated( table_HGNCids$input ) ]
duplicatedinputids_positions <- which(table_HGNCids$input %in% duplicatedinputids)
inputids_that_map_to_more_than_1_HGNC <- table_HGNCids[ duplicatedinputids_positions , ]
table(inputids_that_map_to_more_than_1_HGNC <- table_HGNCids[ duplicatedinputids_positions , ]$input)

# save this table besides the results

# Identifying the HGNCs that received more than one map
duplicatedHGNCsids <- table_HGNCids$name[ duplicated( table_HGNCids$name ) ]
duplicatedHGNCs_positions <- which(table_HGNCids$name %in% duplicatedHGNCsids)
table_of_duplicated_HGNCs <- table_HGNCids[duplicatedHGNCs_positions, ]
table(table_of_duplicated_HGNCs <- table_HGNCids[duplicatedHGNCs_positions, ]$name)
# The distribution of your duplicates = table(duplicatedHGNCsids)

######
## Bulding the matrix with the annotations (HGNC)
######
input_ids_and_converted_ids <- table_HGNCids[ , c(2,4)]

# creating the matrix to fill out
exp_values_part <- matrix( rep(0, dim( input_ids_and_converted_ids)[1] * dim(inputmatrix)[2] ) , ncol = dim(inputmatrix)[2] ) 
converted_ids_uncollapsed_table <- cbind( input_ids_and_converted_ids, exp_values_part) 
colnames(converted_ids_uncollapsed_table)[ 3:dim(converted_ids_uncollapsed_table)[2] ] <- colnames(inputmatrix  )

#######
## filling the matrix section
######
## Take it to another place and make it a funcion
### Maybe using Parallel https://nceas.github.io/oss-lessons/parallel-computing-in-r/parallel-computing-in-r.html
### https://www.r-bloggers.com/2018/09/simple-parallel-processing-in-r/
for( w in 1:dim( converted_ids_uncollapsed_table )[1] ){
  +     position_ENSG_in_inputmat <- grep( converted_ids_uncollapsed_table[w,"input"] , rownames(inputmatrix) )
  +     converted_ids_uncollapsed_table[w , 3:dim(converted_ids_uncollapsed_table)[2] ]   <- inputmatrix[ position_ENSG_in_inputmat, ]
  +   }  
# this can take a long... ~ 1 hr in the server but it was successful

#values_perid <- function(x){
#  position_in_inputmat <- grep(x,rownames(inputmatrix) )
#  inputmatrix[position_in_inputmat,]
#}

# Reading the previous result from anotehr point bc I ran it in the server
#converted_ids_uncollapsed_table_read_readed <-read.table("/media/rmejia/mountme88/Projects/Phosoholipidosis/RNAseq/Expression_Matrix_from_Emmi/lipidosis_RNA_16_STAR_fC_edgeR_matrix_2021_04_23_converted_HGNC_uncollapsed.tsv")

#######
### Collapsing
#######
Collapsed_by_median_dplyr <- converted_ids_uncollapsed_table_read_readed %>%
  group_by(target)%>% 
  summarise(MLipChl1=median(LipChl1), MLipChl2=median(LipChl2), MLipChl3=median(LipChl3), MLipChl4=median(LipChl4),
            MLipCon1=median(LipCon1), MLipCon2=median(LipCon2), MLipCon3=median(LipCon3), MLipCon4=median(LipCon4),
            MNorChl1=median(NorChl1), MNorChl2=median(NorChl2), MNorChl3=median(NorChl3), MNorChl4=median(NorChl4),
            MNorCon1=median(NorChl1), MNorCon2=median(NorChl2), MNorCon3=median(NorChl3), MNorCon4=median(NorChl4)
            )

Collapsed_by_median_dplyr_4_results <- as.data.frame( Collapsed_by_median_dplyr)

# erasing the "target" column and moving it to the "rownames"
rownames(Collapsed_by_median_dplyr_4_results) <- Collapsed_by_median_dplyr_4_results$target
Collapsed_by_median_dplyr_4_results <- Collapsed_by_median_dplyr_4_results[,-grep("target", colnames( Collapsed_by_median_dplyr_4_results))]

###########################
#   Saving the results  ###
###########################
write.table( Collapsed_by_median_dplyr_df, file = path2save , row.names = FALSE, sep="\t" )


