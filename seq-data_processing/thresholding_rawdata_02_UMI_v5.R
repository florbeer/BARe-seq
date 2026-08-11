require(data.table)


#reading in file with metrics of experiment seq run
experiment <- "experiment_name"
sample_path <- paste0("/groups/stark/lorbeer/",experiment,"/Rdata/")
setwd(sample_path)
#file for determining thresholds for each experiment for filtering data from script01
#thresholds <- fread(paste0("/groups/stark/lorbeer/output_directory/out/",experiment,"_raw_qc","/thresholding_",experiment,".txt"))
# determining the threshold for each file based on the cutoff defined in params
# thresholds$cutoff <- round(cutoff * thresholds$median)
# thresholds <- thresholds[which(thresholds$threshold_used ==1),]

#set the output directory
dir <- "/groups/stark/lorbeer/output_directory/"
qc_dir <- paste0(dir,"out/",experiment,"_raw_qc/")

#params
#filter cutoff set to 2%
cutoff <- 0.02

#file for determining thresholds for each experiment for filtering data from script01
thresholds <- fread(paste0(qc_dir,experiment,"_raw_qc","/thresholding_",experiment,".txt"))

# determining the threshold for each file based on the cutoff defined in params
thresholds$cutoff <- round(cutoff * thresholds$median)
thresholds <- thresholds[which(thresholds$threshold_used ==1),]
thresholds 
write.table(thresholds, paste0(sample_path , experiment, "_thresholds_applied.txt"))

# Path where script01 wrote the annotated per-sample files
annot_dir <- paste0(dir,"out/", experiment, "_raw_qc/annotated_samples/")
stopifnot(dir.exists(annot_dir))

#define the sample names based on files in directory 
sample_names <- list.files(annot_dir, pattern = "_annotated\\.tsv\\.gz$", full.names = TRUE)
sample_names

#if it doesn't exist yet, create a directory named
if(!dir.exists(paste0(sample_path,"cutoffv4_",cutoff,"_UMImatch/")))
{
  dir.create(paste0(sample_path,"cutoffv4_",cutoff,"_UMImatch/"))
}
#set as intut directory
out_dir <- paste0(sample_path,"cutoffv4_",cutoff,"_UMImatch/")


#function to filter samples 

sample_cutoff <- function(x){
  # x is a full path to ".../<sample>_annotated.tsv.gz"
  sample <- sub("_annotated\\.tsv\\.gz$", "", basename(x))
  
  # read only this sample (annotated)
  df <- fread(x)
  df <- as.data.table(df)[read_QC %in% c("perfect_read", "mismatch_UMI")]
  
  # get cutoff for this sample from thresholds table
  x_cutoff <- thresholds$cutoff[which(thresholds$sample == sample)]
  if (length(x_cutoff) != 1 || is.na(x_cutoff)) {
    stop(paste0("No unique cutoff found for sample '", sample, "' in thresholds table."))
  }
  
  # count representation per (seq_ID, read_UMI, plasmid_UMI)
  dat_n <- df[,.N, .(seq_ID, read_UMI, plasmid_UMI)]
  
  # keep only high-representation molecules (note your original: N > x_cutoff + 1)
  dat_n <- dat_n[N > (x_cutoff + 1)]
  
  # join N back onto df, then keep rows where N is present
  df2 <- merge(df, dat_n[, .(seq_ID, read_UMI, plasmid_UMI, N)],
               all.x = TRUE, by = c("seq_ID", "read_UMI", "plasmid_UMI"))
  
  # write out using the ORIGINAL raw-like filename (without _annotated...)
  out_file <- file.path(out_dir, paste0(sample, ".txt"))
  
  # your original wrote all columns except the last (which was N after merge)
  # df2 has N as last column from merge, so mimic that behavior:
  keep <- df2[!is.na(N)]
  write.table(keep[, 1:(ncol(keep)-1)],
              out_file, sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

lapply(sample_names, sample_cutoff)
