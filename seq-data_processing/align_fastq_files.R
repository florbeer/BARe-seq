# R version 4.1.2 
require(data.table)
require(Rsubread)
require(stringr)
require(Rsamtools)

#This script is used to align the raw fastq files and create the first version of "counts.txt" files containing all information of the sequencing files 

#set working directory to experiment directory containing subfolder with fastq files 

wdir <- "/groups/stark/lorbeer/experiment_name/"
setwd(wdir)

#read in txt file containing information of indeces for alignment and UMI structure 
# columns "library", "kind", "regel", "pUMI", "pUMI_length", "bowtie_index", "bowtie_index_length", "pattern"    
libfl <- read.table("/groups/stark/lorbeer/analysis_TWIST14_Nov22/data/libFL_table.txt", header = TRUE)

# provide path to custom index for alignment
path_to_index <- "/groups/stark/lorbeer/analysis_TWIST14_Nov22/data/index_twist/"

#------------
# alignment
# first create directories 
if(!dir.exists("./sam/")){
  dir.create("./sam/")
}

if(!dir.exists("./Rdata/")){
  dir.create("./Rdata/")
}

#------------------------------------------------#

#optional step: trimming only required for PE150 reads 
file.edit("/groups/stark/lorbeer/experiment_name/trim_UMI.R")


# identify the files to run 
file.list <-list.files(wdir, recursive = T, full.names = F, pattern = "*.fq.gz")
file_name <- unlist(unique(strsplit(basename(file.list[grepl(pattern , file.list)]), split = "_[1-2].fq.gz")))
#or use trimmed files
lib_file <- sapply( str_split(file_name, pattern = "_"), `[`, 1)
file_table <- data.frame(file_name, lib_file)
file_table <- merge(file_table, libfl, by.x = "lib_file", by.y = "library", all.x = TRUE)

# this file table summarizes the information required for aligning each fastq file of the experiment 
file_table

#---------------
#using the pUMI trimmed files 

alignSTARR <- function(x){
  if(!file.exists(paste0("./sam/", x))){
    Rsubread::align(index = paste0(path_to_index,subset(file_table, file_name==x)$bowtie_index),
          readfile1 = paste0("./fastq/", x, "_1.fq.gz"),
          readfile2 = paste0("./fastq/", x, "_2.fq.gz"),
          type = "dna", 
          input_format = "gzFASTQ", 
          output_file = paste0("./sam/", x), 
          output_format = "SAM",
          maxMismatches = 1,
          #maxMismatches = 3 for 150PE
          unique = T, 
          nTrim5 = 0, 
          nTrim3 = 0, 
          nthreads = 8)
    print(paste0(x, ".sam DONE!"))
  }
}

lapply(file_table$file_name, alignSTARR)


#------------------------------------------------#
#creating a file with the read counts 

#new attempt including mapping stats 
dir.create("./sam_stats", showWarnings = FALSE)

readcounts <- function(x){
  if(!file.exists(paste0("./Rdata/counts_file_",x,".txt"))){
    #read sam file
    # Define paths
    sam_path <- paste0("./sam/", x)
    bam_path <- paste0("./sam/", tools::file_path_sans_ext(x), ".bam")
    
    # Check whether BAM already exists
    if(file.exists(bam_path)){
      message("Using existing BAM file: ", bam_path)
      sam_raw <- scanBam(bam_path)
    } else {
      message("BAM not found, converting SAM → BAM for ", x)
      sam_raw <- scanBam(asBam(sam_path))
    }
    sam <- data.table(read_ID = sam_raw[[1]]$qname, 
                      seq_ID = sam_raw[[1]]$rname, 
                      seq = as.character(sam_raw[[1]]$seq))
    print(paste0(x, " sam read!"))
    
    # Generate SAM stats and save to text file
    generate_sam_stats(x, sam_raw, sam)
    
    #read second read
    fa_2 <- fread(paste0("./", "fastq/", x,"_2.fq.gz"), fill = T, header= F)
    fa_2[, col:= rep(c("read_ID", "rev_read", "V3", "V4"), nrow(fa_2)/4)]
    fa_2[, line:= rep(seq(nrow(fa_2)/4), each= 4)]
    fa_2 <- dcast(fa_2, line~col, value.var = "V1")
    #make read_ID match with sam file
    fa_2[, read_ID:= gsub("^@","", read_ID)]
    print(paste0(x, " fa_2 read!"))
    
    res <- merge(sam, 
                 fa_2[, .(read_ID, rev_read)], 
                 by= "read_ID")
    res[, plasmid_UMI:= substr(rev_read, 1, subset(file_table, file_name==x)$pUMI_length)]
    res[, read_UMI:= tstrsplit(read_ID, "_", keep= 2)]
    res$read_ID <- NULL
    res[,handle := substr(rev_read, subset(file_table, file_name==x)$pUMI_length+1,36)]
    fwrite(res, paste0("./Rdata/counts_file_",x,".txt"), col.names= T, row.names= F, sep= "\t", quote= F)
    print(paste0(x, " DONE!"))
  }
}

# Helper function to generate SAM stats
generate_sam_stats <- function(sample_name, sam_raw, sam_dt) {
  # Calculate basic alignment stats
  total_reads <- length(sam_raw[[1]]$qname)
  mapped_reads <- sum(!is.na(sam_raw[[1]]$rname))
  unmapped_reads <- total_reads - mapped_reads
  alignment_rate <- round(100 * mapped_reads / total_reads, 2)
  
  # Multi-mapped reads (approx. via MAPQ=0)
  multi_mapped_reads <- sum(sam_raw[[1]]$mapq == 0, na.rm = TRUE)
  unique_mapped_reads <- mapped_reads - multi_mapped_reads
  multi_rate <- round(100 * multi_mapped_reads / total_reads, 2)
  
  # Count reads per reference sequence
  ref_counts <- table(sam_raw[[1]]$rname, useNA = "ifany")
  ref_counts <- ref_counts[order(ref_counts, decreasing = TRUE)]
  
  # Calculate mapping quality stats
  mapq_stats <- summary(sam_raw[[1]]$mapq[!is.na(sam_raw[[1]]$mapq)])
  
  # Create stats summary
  stats_text <- c(
    paste0("SAM Alignment Statistics for: ", sample_name),
    paste0("Generated on: ", Sys.time()),
    "",
    "=== BASIC STATS ===",
    paste0("Total reads: ", total_reads),
    paste0("Mapped reads: ", mapped_reads),
    paste0("  • Unique-mapped reads: ", unique_mapped_reads),
    paste0("  • Multi-mapped reads: ", multi_mapped_reads,
           " (", multi_rate, "% of total)"),
    paste0("Unmapped reads: ", unmapped_reads),
    paste0("Alignment rate: ", alignment_rate, "%"),
    "",
    "=== MAPPING QUALITY STATS ===",
    paste0("Min MAPQ: ", mapq_stats["Min."]),
    paste0("Mean MAPQ: ", round(mapq_stats["Mean"], 2)),
    paste0("Max MAPQ: ", mapq_stats["Max."]),
    "",
    "=== READS PER REFERENCE ===",
    paste0(names(ref_counts), ": ", ref_counts, collapse = "\n")
  )
  
  # Write stats to file
  writeLines(stats_text, paste0("./sam_stats/", sample_name, "_alignment_stats.txt"))
  print(paste0("SAM stats saved for ", sample_name))
}

lapply(file_table$file_name, readcounts)

