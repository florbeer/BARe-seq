require(vroom)
require(data.table)
require(ggplot2)
library(moments)
library(dplyr)
require(stringr)
require(patchwork)
require(ggpubr)
require(purrr)

#reading in files
experiment <- "experiment_name"
sample_path <- paste0("/groups/stark/lorbeer/",experiment,"/Rdata/")
setwd(sample_path)
#cutoff <- 0.02
sample_names <- list.files(sample_path, pattern = "counts_file_lib.*\\.txt$")

#read in txt file containing information of indeces for alignment and UMI structure 
# columns "library", "kind", "regel", "pUMI", "pUMI_length", "bowtie_index", "bowtie_index_length", "pattern"    
libfl <- read.table("/groups/stark/lorbeer/analysis_TWIST14_Nov22/data/libFL_table.txt", header = TRUE)

existing <- file.exists(paste0(out_dir, sample_names))
if (any(!existing)) {
  message("Dropping missing combined files: ", paste(sample_names[!existing], collapse = ", "))
}
sample_names <- sample_names[existing]

#directory where collapsed files can be found
in_dir <- paste0("/groups/stark/lorbeer/",experiment,"/Rdata/cutoffv4_0.02_UMImatch/log_enhancer_submission/lv_tUMI_collapse/")

#if it doesn't exist yet, create a directory named UMIcollapsed
if(!dir.exists(paste0("//groups/stark/lorbeer/output_directory/out/",experiment,"_processed/")))
{
  dir.create(paste0("/groups/stark/lorbeer/output_directory/out/",experiment,"_processed/"))
}

#set as iytut directory
out_dir <- paste0("/groups/stark/lorbeer/output_directory/out/",experiment,"_processed/")


#This function combines the files and saves them in the out_dir one file per sample
combine_aftercollapse <- function(x, in_dir, out_dir){
  
  output_file <- paste0(out_dir, x)
  # Check if the output file already exists to avoid overwriting
  if (!file.exists(output_file)) {
    files <- list.files(in_dir, pattern = toString(paste0(sub('\\.txt$', '', x),".*lv3collapsed.txt_main_table.txt_tUMI_collapsed.txt$")), full.names = "FALSE")
    dat_collapsed <- rbindlist(lapply(paste0(in_dir,files), fread))
    dat_collapsed$RNA_molecule <- paste0(dat_collapsed$seq_ID, "_", dat_collapsed$plasmid_UMI_main, "_", dat_collapsed$read_UMI_main)
    fwrite(dat_collapsed, output_file, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    message("File saved: ", output_file)
  } else {
    message("File already exists, skipping: ", output_file)
  }
}

lapply(sample_names, combine_aftercollapse, in_dir, out_dir)

#information on libraries used with libfl file
file_name <- c()
for (i in (1:length(sample_names))){
  x <- sample_names[i]
  sample <-basename(substr(toString(x),13, nchar(x)-4))
  file_name[i]<- sample
}
lib_file <- sapply( str_split(file_name, pattern = "_"), `[`, 1)
file_table <- data.frame(file_name, lib_file)
file_table <- merge(file_table, libfl, by.x = "lib_file", by.y = "library", all.x = TRUE)

# where to store per-sample intermediate outputs
counts_dir <- paste0(out_dir, "per_sample_counts/")
qc_dir <- paste0(out_dir, "per_sample_qc/")
if (!dir.exists(counts_dir)) dir.create(counts_dir, recursive = TRUE)
if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

# safer summary function (fixes your dat2 column name mismatch)
summarize_qc_dt <- function(dt, name) {
  counts <- dt[, .N, .(sample, seq_ID, plasmid_UMI_main, UMI_match)]
  dat1 <- counts[, .N, .(sample, seq_ID)]
  dat2 <- counts[, .(No_transcripts = sum(N)), .(sample, seq_ID)]
  dat3 <- merge(dat1, dat2, by = c("sample", "seq_ID"))
  dat3[, Normalized_counts_sequence := No_transcripts / N]
  
  data.frame(
    Name = name,
    total_unique_RNA_molecules = nrow(unique(dt)),
    total_unique_plasmids = length(unique(dt$plasmid_UMI_main)),
    Median_Amplification_counts_final = median(dt$Amplification_counts_final, na.rm = TRUE),
    Mean_Amplification_counts_final = mean(dt$Amplification_counts_final, na.rm = TRUE),
    median_pUMI_perseq = median(dat1$N),
    median_tx_perseq = median(dat2$No_transcripts),
    median_normalized_txperplasmid = median(dat3$Normalized_counts_sequence)
  )
}

process_one_sample <- function(x, out_dir, file_table, counts_dir, qc_dir) {
  # x is something like "counts_file_lib...txt"
  in_file <- paste0(out_dir, x)
  if (!file.exists(in_file)) stop(paste0("Missing combined file: ", in_file))
  
  dt <- as.data.table(fread(in_file))
  
  sample <- basename(substr(toString(x), 13, nchar(x) - 4))
  
  pat <- file_table[which(file_table$file_name == sample), ]$pattern
  if (length(pat) != 1 || is.na(pat)) stop(paste0("No unique pattern for sample: ", sample))
  
  dt[, sample := sample]
  dt[, UMI_match := grepl(pat, plasmid_UMI_main)]
  
  # 1) per-sample counts table saved to disk (used later for global rbind)
  counts <- dt[, .N, .(sample, seq_ID, plasmid_UMI_main, UMI_match)]
  fwrite(counts,
         file = paste0(counts_dir, x, ".counts.tsv.gz"),
         sep = "\t")
  
  # 2) per-sample QC summary row saved to disk
  qc_row <- summarize_qc_dt(dt, x)
  fwrite(as.data.table(qc_row),
         file = paste0(qc_dir, x, ".qc.tsv"),
         sep = "\t")
  
  # 3) return the dt ONLY for immediate plotting (not stored globally)
  dt
}

#----
#QC table 
#-----

pdf(paste0(out_dir, experiment, "_final_Ampcounts_histograms.pdf"))

for (x in sample_names) {
  in_file <- paste0(out_dir, x)
  
  if (!file.exists(in_file)) {
    message("SKIP (missing combined file): ", in_file)
    next
  }
  
  message("Plotting: ", x)
  
  dt <- process_one_sample(x, out_dir, file_table, counts_dir, qc_dir)
  
  dt[, Amplification_counts_final := as.numeric(Amplification_counts_final)]
  
  p <- ggplot(dt) +
    geom_histogram(aes(x = Amplification_counts_final, fill = UMI_match),
                   binwidth = 1, na.rm = TRUE) +
    theme_classic() +
    ggtitle(x)
  
  print(p)
  
  rm(dt)
  gc()
}

dev.off()

# qc part 
qc_files <- list.files(qc_dir, pattern = "\\.qc\\.tsv$", full.names = TRUE)
qc_summary <- rbindlist(lapply(qc_files, fread), fill = TRUE)

write.table(qc_summary,
            file = paste0(out_dir, experiment, "_qc_summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(qc_summary)

#now create a table with all criteria summarized
count_files <- list.files(counts_dir, pattern = "\\.counts\\.tsv\\.gz$", full.names = TRUE)
counts <- rbindlist(lapply(count_files, fread), fill = TRUE)


#----
#next section 
#---- 

#now create a table with all criteria summarized
counts <- rbindlist(dat)[, .N, .(sample, seq_ID, plasmid_UMI_main, UMI_match)]

#saveRDS(counts, paste0(out_dir,experiment,"all_counts.rds"))
#counts <- readRDS(paste0(out_dir,experiment,"_counts.rds"))

counts <- counts %>% filter(UMI_match == TRUE)

# save the table used for downstream analyses and provided as processed data table at GEO 
write.table(counts, paste0(out_dir,experiment,"_processing1_counts.txt"), sep="\t", row.names = F, quote=F)

#---------------
#barcode QC plots according to Reyna 
#Plot No 1: Number of tUMIs/plasmidUMI (Reyna N per plasmid_UMI_main)

#barcodes per sequence 
dat1 <- counts[, .N, .(sample, seq_ID, experiment, library, extraction, replicate)]
p <- ggplot(dat1) + 
  geom_boxplot(aes(y = N, x = sample, color = extraction, fill = replicate)) + 
  #scale_fill_manual(values = c("lightblue", "darkblue")) + 
  #scale_color_manual(values = c("darkred", "coral2")) + 
  scale_y_log10() + 
  #coord_cartesian(ylim=c(0.99,2.5)) + 
  scale_x_discrete(guide = guide_axis(angle =45)) +
  ggtitle("barcodes per regulatory element") + 
  #facet_wrap(~library, ncol = 1) + 
  theme_classic()

#counts per sequence 
dat2 <- counts[, .(No_transcripts = sum(N)), .(sample, seq_ID, experiment, library, extraction, replicate)]
q <- ggplot(dat2) + 
  geom_boxplot(aes(y = No_transcripts, x = sample, color = extraction, fill = replicate)) + 
  #scale_fill_manual(values = c("lightblue", "darkblue")) + 
  #scale_color_manual(values = c("darkred", "coral2")) + 
  scale_y_log10() + 
  #coord_cartesian(ylim=c(0.99,2.5)) + 
  scale_x_discrete(guide = guide_axis(angle =45)) +
  ggtitle("transcript counts per regulatory element") + 
  #facet_wrap(~library, ncol = 1) + 
  theme_classic()

#normalized counts per sequence 
dat3 <- merge(dat1, dat2)
dat3[,Normalized_counts_sequence := No_transcripts/N, .(sample, seq_ID, extraction, replicate)]
r <- ggplot(dat3) + 
  geom_boxplot(aes(y = Normalized_counts_sequence, x = sample, color = extraction, fill = replicate)) + 
  #scale_fill_manual(values = c("lightblue", "darkblue")) + 
  #scale_color_manual(values = c("darkred", "coral2")) +
  scale_y_log10() + 
  #coord_cartesian(ylim=c(0.99,2.5)) + 
  scale_x_discrete(guide = guide_axis(angle =45)) +
  ggtitle("normalized counts per sequence (transcript/plasmid per reg element)") + 
  #facet_wrap(~library, ncol = 1) + 
  theme_classic()

figure <- ggarrange(p,q,r, ncol = 1)
ggexport(figure, filename = paste0(out_dir, experiment, "_qc_boxplots.pdf"), height = 16)



