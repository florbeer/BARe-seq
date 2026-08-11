options(stringsAsFactors=FALSE)
options("scipen"=100, "digits"=4)

#test options
opt=list()
opt$input <- "/groups/stark/lorbeer/E21_TWIST14_STAP-seq/Rdata/cutoff_0.02/log_enhancer_submission/counts_file_libFL008_E21_kf_1_m2_chr2R_9812941_9813074_-_lv3collapsed.txt_main_table.txt"
opt$MM <- 1
opt$core <- 4
opt$out <- "counts_file_libFL008_E21_kf_1_m2_chr2R_9812941_9813074_-_lv3collapsed.tUMI.txt"

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-i", "--input"),  action = "store", type="character", default=NULL,
              help="filtered_libFL001_E04_S1-2_control_flat_genomic_region_A_00702_collapsed.txt_main_table.txt", metavar="character"),
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
  description = "tUMI collapsing"
)
arguments <- parse_args(opt_parser, positional_arguments = TRUE)
opt <- arguments$options

#------------
# LIBRARIES
#------------

require(data.table)
require(dplyr)
require(stringr)
require(stringdist)
require(parallel)

#######

#------------
# Prepare data
#------------

test_big <- fread(opt$input)
test_big$enh_pUMI_tUMI_ID <- paste(test_big$seq_ID, test_big$plasmid_UMI_main, test_big$read_UMI, sep="_")
t <- data.frame(table(test_big$enh_pUMI_tUMI_ID))
test_big$moleculeCounts <- test_big$Amplification_counts + t$Freq[match(test_big$enh_pUMI_tUMI_ID, t$Var1)] - 1 # remove its own counting

test_big_sorted <- setorder(test_big, plasmid_UMI_main, -moleculeCounts)

f3_tUMI = function(bar,c=1){
  keep <- bar
  while (length(bar)>0) {
    rmv = which(stringdist(bar[1],bar,method="hamming",nthread =c)<=opt$MM)
    keep[names(bar)[rmv]] = bar[1] # replace the main tUMI in the place of collapsed ones
    bar = bar[-rmv]
  }
  return(keep)
}


#------------
# Run
#------------

# Calculate the number of cores
no_cores <- opt$core
# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl, "stringdist")
clusterExport(cl, "f3_tUMI")
clusterExport(cl, "opt")

results <- unlist(parLapply(cl, split(test_big_sorted$read_UMI,test_big_sorted$plasmid_UMI_main),
                            function(x){
                              names(x) <- 1:length(x)
                              f3_tUMI(x)
                            })
)

# Finish
stopCluster(cl)

test_big_sorted$read_UMI_main <- as.character(results)
# print to see top
head(test_big_sorted)
#test_big_sorted <- test_big_sorted[test_big_sorted$counts>0]

# get final amplification counts
test_big_sorted$enh_pUMI_tUMI_ID <- paste(test_big_sorted$seq_ID, test_big_sorted$plasmid_UMI_main, test_big_sorted$read_UMI_main, sep="_")
t <- data.frame(table(test_big_sorted$enh_pUMI_tUMI_ID))
test_big_sorted$Amplification_counts_2 <- t$Freq[match(test_big_sorted$enh_pUMI_tUMI_ID, t$Var1)]

test_big_sorted$Amplification_counts_final <- test_big_sorted$Amplification_counts + test_big_sorted$Amplification_counts_2 -1

# remove duplicated molecules
out <- test_big_sorted[!duplicated(test_big_sorted$enh_pUMI_tUMI_ID),]

# print output
write.table(out[,c(2,3,1,4,5,7,8,9,10)], opt$out, sep="\t", row.names = F, quote = F)


#### plots
# t <- data.frame(table(out$plasmid_UMI_main))
# out$RNAcounts <- t$Freq[match(out$plasmid_UMI_main, t$Var1)]


