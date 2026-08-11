#!/usr/bin/env Rscript

# setwd("/groups/stark/almeida/Projects/Franzi_bursting/results/20200131_test_UMI_collapsing_script/")

# run as
# cd /groups/stark/almeida/Projects/Franzi_bursting/results/20200131_test_UMI_collapsing_script/
# Rexec="singularity run --app Rscript /groups/stark/software-all/singularity_img/singularity.R.with_pckg.simg "
# bsub "$Rexec /groups/stark/almeida/Projects/Franzi_bursting/src/Franzi_pUMI_collapsing_submission_jobs.R -i /groups/stark/lorbeer/E01_STARR_burst_reseq/Rdata/filtered_table_libFL001_STARRseq_01.txt -t log_enhancer_submission"

options(stringsAsFactors=FALSE)
options("scipen"=100, "digits"=4)

#test options
opt=list()
opt$input <- "/groups/stark/lorbeer/E01_STARR_burst_reseq/Rdata/filtered_table_libFL001_STARRseq_01.txt"
opt$tmp_folder <- "log_enhancer_submission"

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-i", "--input"),  action = "store", type="character", default=NULL,
              help="final_table_libFL001_STARRseq_S1.1", metavar="character"),
  make_option(c("-t", "--tmp_folder"),  action = "store", type="character", default="log_enhancer_submission",
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

opt

#------------
# Load data
#------------

cat("\nReading main table\n")
test_big <- read.delim(opt$input, stringsAsFactors = F)

#------------
# submit job
#------------

# create tmp directory
if(!dir.exists(opt$tmp_folder)) system(paste0("mkdir ", opt$tmp_folder))

prefix <- substr(basename(opt$input), 1, nchar(basename(opt$input))-4)

# loop per enhancer
for(enh in unique(test_big$seq_ID)){
  
  # subset table
  tmp <- test_big[test_big$seq_ID %in% enh,]
  
  # write reads table
  write.table(tmp, paste0(opt$tmp_folder, "/", prefix, "_", enh, ".txt"), sep="\t", row.names = F, quote = F)
  
  # define job time
  if(nrow(tmp)>50000){
    t="1-00:00:00"
  }else{
    t="1:00:00"
  }
  # submit job
  system(paste0("/groups/stark/software-all/shell/bsub_gridengine -o ", opt$tmp_folder, " -T ", t,
                " 'singularity run --app Rscript /groups/stark/software-all/singularity_img/singularity.R.with_pckg.simg /groups/stark/lorbeer/functions/Franzi_pUMI_collapsing_LSdistance_amp_counts_v2.R -i ", opt$tmp_folder, "/", prefix, "_", enh, ".txt -m 3 -c 1 -o ", opt$tmp_folder, "/", prefix, "_", enh, "_lv3collapsed.txt' > ", opt$tmp_folder, "/msg.", enh, ".tmp"))
  
  print(paste0("summited ", enh, " job"))
}

sessionInfo()
