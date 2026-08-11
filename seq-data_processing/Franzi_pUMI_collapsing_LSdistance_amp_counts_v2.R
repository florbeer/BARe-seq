#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
options("scipen"=100, "digits"=4)

#test options
opt=list()
opt$input <- "/groups/stark/lorbeer/E21_TWIST14_STAP-seq/Rdata/cutoff_0.02/counts_file_libFL002_E21_kf_1_m2.txt"
opt$MM <- 3
opt$core <- 4
opt$out <- "counts_file_libFL002_E21_kf_1_m2_collapsed.txt"

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-i", "--input"),  action = "store", type="character", default=NULL,
              help="final_table_libFL001_STARRseq_S1.1", metavar="character"),
  make_option(c("-m", "--MM"),  action = "store", type="numeric", default=1,
              help="pUMI_MM", metavar="character"),
  make_option(c("-c", "--core"),  action = "store", type="numeric", default=1,
              help="number of cores", metavar="character"),
  make_option(c("-o", "--out"),  action = "store", type="character", default="out.txt",
              help="out txt file", metavar="character")
)

opt_parser <- OptionParser(
  usage = "%prog [options]",
  option_list=option_list,
  description = "pUMI collapsing"
)
arguments <- parse_args(opt_parser, positional_arguments = TRUE)
opt <- arguments$options

#------------
# LIBRARIES
#------------

suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("stringdist"))

# print options
cat("\nRunning UMI collapsing\n")

opt

#------------
# Prepare data
#------------

test_big <- read.delim(opt$input, stringsAsFactors = F)
test_big$ID <- paste(test_big$seq_ID, test_big$plasmid_UMI, test_big$read_UMI, sep="_")
t <- data.frame(table(test_big$ID))
test_big$Amplification_counts <- t$Freq[match(test_big$ID, t$Var1)]

# collapse by enhancer-pUMI-rUMI (uniquely)
test_big <- test_big[!duplicated(test_big$ID),]

# get RNA counts per enhancer-pUMI just to order them and find the main one
test_big$enh_pUMI_ID <- paste(test_big$seq_ID, test_big$plasmid_UMI, sep="_")
t <- data.frame(table(test_big$enh_pUMI_ID))
test_big$RNA_counts <- t$Freq[match(test_big$enh_pUMI_ID, t$Var1)]

# function to get main pUMI
f3_pUMI = function(bar,c=1){
  keep <- bar
  while (length(bar)>0) {
    rmv = which(stringdist(bar[1],bar,method="lv",nthread =c)<=opt$MM)
    keep[names(bar)[rmv]] = bar[1] # replace the main pUMI in the place of collapsed ones
    bar = bar[-rmv]
  }
  return(keep)
}

#------------
# replace the main pUMI in the place of the collapsed ones
#------------

cat("\nReplacing main pUMI at coolapsed ones\n")

# Calculate the number of cores
no_cores <- opt$core
# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl, "stringdist")
clusterExport(cl, "f3_pUMI")
clusterExport(cl, "opt")

tmp <- test_big
tmp <- tmp[order(tmp$RNA_counts, decreasing = T),] # order by counts to get the main pUMI on top
# tmp <- test_big[test_big$seq_ID %in% unique(test_big$seq_ID)[1:5],]
tmp <- tmp[order(tmp$seq_ID),] # order by enhancer

### Get collpased IDs
clusterExport(cl, "tmp")
results <- unlist(parLapply(cl, split(tmp$plasmid_UMI,tmp$seq_ID),
                                        function(x){
                                          names(x) <- 1:length(x)
                                          f3_pUMI(x)
                                        })
)


tmp$plasmid_UMI_main <- as.character(results)

#### Maybe we could just output here the table with enhancer, pUMI (new, after collapsing) and rUMI - because then this table is like the initial one
write.table(tmp[,c("seq_ID", "read_UMI","plasmid_UMI", "plasmid_UMI_main", "Amplification_counts")], paste0(opt$out, "_main_table.txt"), sep="\t", row.names = F, quote = F)

sessionInfo()
