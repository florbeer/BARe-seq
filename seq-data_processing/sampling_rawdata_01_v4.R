#load packages
require(vroom)
require(data.table)
require(stringr)

#-----------------
# parameters of the analysis 
#----------------

#change the experiment name 
experiment <- "experiment_name"
#set the output directory
dir <- "/groups/stark/lorbeer/output_directory/"

#sample_path to output folder of previous script
sample_path <- paste0("/groups/stark/lorbeer/",experiment,"/Rdata/")
setwd(sample_path)

#directory where the results of this analysis are saved
if(!dir.exists(paste0(dir,"out/",experiment, "_raw_qc")))
{
  dir.create(paste0(dir,"out/",experiment,"_raw_qc"))
}
out_dir <- paste0(dir,"out/",experiment,"_raw_qc/")

#read in txt file containing information of indeces for alignment and UMI structure 
# columns "library", "kind", "regel", "pUMI", "pUMI_length", "bowtie_index", "bowtie_index_length", "pattern"    
libfl <- read.table("/groups/stark/lorbeer/analysis_TWIST14_Nov22/data/libFL_table.txt", header = TRUE)
#and list the files that will be analyzed for this experiment
sample_names <- list.files(sample_path, pattern = "counts_file_lib.*rep[1-7]\\.txt$")
sample_names

#information on libraries used is combined for each file
file_name <- c()
for (i in (1:length(sample_names))){
  x <- sample_names[i]
  sample <-basename(substr(toString(x),13, nchar(x)-4))
  file_name[i]<- sample
}
lib_file <- sapply( str_split(file_name, pattern = "_"), `[`, 1)
file_table <- data.frame(file_name, lib_file)
file_table <- merge(file_table, libfl, by.x = "lib_file", by.y = "library", all.x = TRUE)

#------------
#annotate file 
#------------
annotate_one_sample <- function(infile, file_table, out_annot_dir) {
  # infile is the filename from sample_names (e.g. "counts_file_lib...txt")
  
  # derive sample key exactly like you do elsewhere
  name <- basename(substr(toString(infile), 13, nchar(infile) - 4))
  
  # pull the pattern for this sample
  pat <- file_table[which(file_table$file_name == name), ]$pattern
  if (length(pat) != 1 || is.na(pat)) {
    stop(paste0("No unique pattern found for sample '", name, "' (check file_table merge)."))
  }
  
  # read only this file
  x <- vroom(infile, col_types = c("c","c","c","c","c","c"))
  x <- as.data.table(x)
  
  # compute QC columns (same logic as your original)
  x[, UMI_match := grepl(pat, plasmid_UMI) & !grepl("GGGGGG", plasmid_UMI)]
  x[, not_aligned := grepl("\\*", seq_ID) | is.na(seq_ID)]
  x[, read_QC := fifelse(UMI_match & !not_aligned, "perfect_read",
                         fifelse(!UMI_match & !not_aligned, "mismatch_UMI",
                                 fifelse(UMI_match &  not_aligned, "not_aligned",
                                         "no_match_atall")))]
  
  # write annotated file for downstream use
  if (!dir.exists(out_annot_dir)) dir.create(out_annot_dir, recursive = TRUE)
  
  out_file <- file.path(out_annot_dir, paste0(name, "_annotated.tsv.gz"))
  fwrite(x, out_file, sep = "\t")
  
  out_file
}


#----
#run annotation for all samples 
# ---- 

out_annot_dir <- file.path(out_dir, "annotated_samples")

annot_files <- vapply(
  sample_names,
  annotate_one_sample,
  FUN.VALUE = character(1),
  file_table = file_table,
  out_annot_dir = out_annot_dir
)

annot_files


#------------------
#make the plots per sample showing the read distribution 
#------------------
require(data.table)
require(ggplot2)
require(gridExtra)

qc_samples_from_file <- function(infile_annot) {
  # infile_annot is like ".../annotated_samples/<sample>_annotated.tsv.gz"
  sample <- sub("_annotated\\.tsv\\.gz$", "", basename(infile_annot))
  
  df <- fread(infile_annot)
  df <- as.data.table(df)
  
  dat_n <- df[,.N, .(seq_ID, read_UMI, plasmid_UMI, UMI_match, not_aligned, read_QC)]
  
  dat_sum <- table(dat_n$N, dat_n$read_QC)
  dat_sum <- cbind(No_reads = c(as.numeric(rownames(dat_sum))), dat_sum)
  colnames(dat_sum) <- c("No_reads", "mismatch_UMI", "no_match_atall", "not_aligned", "perfect_read")
  dat_sum <- as.data.frame(dat_sum)
  
  dat_sum$reads_mismatch_UMI   <- dat_sum$No_reads * dat_sum$mismatch_UMI
  dat_sum$reads_no_match_atall <- dat_sum$No_reads * dat_sum$no_match_atall
  dat_sum$reads_not_aligned    <- dat_sum$No_reads * dat_sum$not_aligned
  dat_sum$reads_perfect_read   <- dat_sum$No_reads * dat_sum$perfect_read
  
  dat_sum_melt <- melt(as.data.table(dat_sum), id.vars = c("No_reads"),
                       measure.vars = c("reads_mismatch_UMI", "reads_no_match_atall",
                                        "reads_not_aligned", "reads_perfect_read"))
  dat_sum_melt[dat_sum_melt == 0] <- NA
  
  reads_correct <- vector("list", 80)
  for(i in 1:80){
    reads_correct[[i]] <- sum(dat_sum[dat_sum$No_reads >= i, ]$mismatch_UMI) /
      (sum(dat_sum[dat_sum$No_reads >= i, ]$perfect_read) +
         sum(dat_sum[dat_sum$No_reads >= i, ]$mismatch_UMI))
  }
  
  cutoff_table <- cbind(
    cutoff = seq_along(reads_correct),
    fraction = unlist(reads_correct)
  )
  write.table(cutoff_table,
              paste0(out_dir, "cutoff_FALSEUMI_", experiment, "_", sample, ".txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  png(file = paste0(out_dir, sample, "_read_plots_UMI.png"), type = "cairo", res = 300,
      width = 9, height = 9, units = "in")
  
  p <- ggplot(dat_n) +
    geom_histogram(aes(x = N, fill = read_QC), binwidth = 1) +
    xlab("N Read representation") +
    scale_y_log10() +
    xlim(0,500) +
    scale_fill_manual(values = c("snow2", "snow3", "snow4", "navyblue")) +
    ggtitle(paste0(sample, " :All raw reads based on sequence matches UMI/RE")) +
    theme_classic()
  
  q <- ggplot(dat_sum_melt) +
    geom_point(aes(x = No_reads, y = value, color = variable)) +
    scale_color_manual(values = c("snow2", "snow3", "snow4", "navyblue")) +
    facet_wrap(~variable) +
    xlab("N Read representation") +
    ylab("Count") +
    xlim(0,500) +
    ggtitle(paste0(sample, " :Sequencing depth of reads by QC")) +
    theme_classic()
  
  r <- ggplot(as.data.frame(cutoff_table)) +
    geom_point(aes(x = cutoff, y = fraction)) +
    geom_hline(yintercept = 0.1, linetype = "dashed") +
    xlab("Cutoff >= x") +
    ylab("Fraction of correct reads") +
    xlim(0,500) +
    ggtitle(paste0(sample, " :Fraction of correct (UMImatch) reads based on representation cutoff")) +
    theme_classic()
  
  grid.arrange(p, q, r, ncol = 1,  layout_matrix = cbind(c(2,2,1), c(2,2,3)))
  dev.off()
  
  print(paste0("qc graphs for ", sample, " have been made!"))
  invisible(cutoff_table)
}

#apply the function
lapply(annot_files, qc_samples_from_file)


#---------------
# making a table for thresholding the samples
#-------------

#thresholds to be applied for characterizing the raw data
vec_threshold <- c(1, 0.75, 0.5, 0.2, 0.1, 0.05, 0.01, 0.005)


sample_sample_from_file <- function(infile_annot) {
  sample <- sub("_annotated\\.tsv\\.gz$", "", basename(infile_annot))
  
  df <- fread(infile_annot)
  df <- as.data.table(df)[read_QC == "perfect_read"]
  
  threshold_info <- list()
  for (n in 1:length(vec_threshold)) {
    threshold <- vec_threshold[n]
    
    # keep your sampling behavior, but guard against tiny files
    take_n <- max(1, floor(threshold * nrow(df)))
    dat_n <- df[sample(nrow(df), take_n), ]
    
    dat_n <- dat_n[,.N, .(seq_ID, read_UMI, plasmid_UMI)]
    dat_n <- as.data.table(cbind(as.numeric(names(table(dat_n$N))), as.vector(table(dat_n$N))))
    colnames(dat_n) <- c("No_reads", "Count")
    dat_n$reads <- dat_n$No_reads * dat_n$Count
    
    threshold_info[[n]] <- data.table(
      sample = sample,
      No_raw_reads = nrow(df),
      No_unique_reads = dim(unique(df[,-1]))[1],
      No_samp_reads = take_n,
      reads_retained = sum(dat_n[No_reads > 1]$Count) / sum(dat_n[No_reads >= 1]$Count),
      threshold_used = threshold,
      median = quantile(dat_n$No_reads, c(.5)),
      high95_percentile = quantile(dat_n$No_reads, c(.95))
    )
  }
  
  rbindlist(threshold_info)
}

info <- rbindlist(lapply(annot_files, sample_sample_from_file))

# save a table with the metrics (will be used to determine the filtering params downstream)
write.table(info, paste0(out_dir,"thresholding_",experiment,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#------------------
#Create the plots showing sequencing depth and read distribution

png(file = paste0(out_dir,experiment,"_threshold_v_readsretained.png"), type = "cairo", res = 300, 
    width = 15, height = 12, units = "in")
par(mfcol = c(2,1))
ggplot(data = info) + 
  geom_line(aes(x = -threshold_used, y = reads_retained)) + 
  geom_point(aes(x = -threshold_used, y = reads_retained)) + 
  facet_wrap(~sample)+
  theme_classic()
dev.off()

png(file = paste0(out_dir,experiment,"_median_v_readsretained.png"), type = "cairo", res = 300, 
    width = 15, height = 12, units = "in")
par(mfcol = c(2,1))
ggplot(data = info) + 
  geom_line(aes(x = log(median), y = reads_retained)) + 
  geom_point(aes(x = log(median), y = reads_retained)) + 
  geom_vline(xintercept = 5, color = "red") + 
  facet_wrap(~sample)+
  theme_classic()
dev.off()

png(file = paste0(out_dir,experiment,"_threshold_v_median.png"), type = "cairo", res = 300, 
    width = 15, height = 12, units = "in")
par(mfcol = c(2,1))
ggplot(data = info) + 
  geom_line(aes(y = threshold_used, x = -log(median))) + 
  geom_point(aes(y = threshold_used, x = -log(median))) + 
  geom_vline(xintercept = -5, color = "red") + 
  facet_wrap(~sample)+
  theme_classic()
dev.off()



