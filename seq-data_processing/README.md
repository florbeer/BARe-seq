# Scripts for processing raw sequencing data (fastq files) to count files required for burst inference scripts
 ## Run in the following order for fastq processing: 
 1) align_fastq_files.R
    - this script aligns and creates *_counts.txt files 
 2) sampling_rawdata_01_v4.R
    - this script generates QC plots for library sequencing depth assessment
    - thresholds.txt file 
 3) thresholding_rawdata_02_UMI_v5.R
    - uses cutoff of 2% of median library coverage to filter reads to accurate ensure counting of real transcripts
    - saves new read coverage filtered counts txt file 
 4) pUMI collapsing: wrapper script: Franzi_pUMI_collapsing_submission_jobs_5lv_amp_v2.R calling collapsing script: Franzi_pUMI_collapsing_LSdistance_amp_counts_v2.R
    - pUMI collapsing, submitted from command line: Within each oligo (seq_ID), lower-abundance pUMIs within a Levenshtein distance of 3 from a higher-abundance pUMI were collapsed utilizing the stringdist package (version 0.9.8).
    - saves one txt file per seq_ID with new columns: "seq_ID"	"read_UMI"	"plasmid_UMI"	"plasmid_UMI_main"	"Amplification_counts"
 5) tUMI collapsing: Franzi_tUMI_collapse_amp_counts.R
    - execute from command line to obtain one txt file per seq_ID with the addditional columns: "moleculeCounts"	"read_UMI_main"	"Amplification_counts_2"	"Amplification_counts_final"
 6) combining_collapsed_reads_05_v3_QC.R
    - combines the files from individual seq_IDs used for downstream analyses (burst inference)
    - produces QC plots to assess the experiment
