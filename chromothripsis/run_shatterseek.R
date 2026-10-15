#!/usr/bin/env Rscript

suppressMessages(library(ShatterSeek))
suppressMessages(library(GenomicRanges))

# =============================================================================
# Usage:
#   Rscript run_shatterseek.R <sample_list_file> [output_name]
#
# Arguments:
#   sample_list_file : tab-separated file with three columns: sampleID, SV_file, CNV_file
#   output_basename  : name for output files (default: "test_shatterseek")
#                      Results will be written as <name>_result.txt and <name>_result_fdr.txt
# =============================================================================

#' @param sv_file Path to structural variant file (tab separated, must contain chrom1, pos1, chrom2, pos2, SVtype, strand1, strand2)
#' @param cnv_file Path to copy number variant file (tab separated, must contain chromosome, start, end, total_cn)
#' @param out_prefix Output file prefix. 

run_shatterseek_sample <- function(sv_file, cnv_file, out_prefix) {
  
  # Read SV data
  sv <- read.table(sv_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  SV_data <- SVs(chrom1 = as.character(sv$chrom1),
                 pos1  = as.numeric(sv$pos1),
                 chrom2 = as.character(sv$chrom2),
                 pos2  = as.numeric(sv$pos2),
                 SVtype = as.character(sv$SVtype),
                 strand1 = as.character(sv$strand1),
                 strand2 = as.character(sv$strand2))
  
  # Read CNV data and merge adjacent segments with identical copy number
  d <- read.table(cnv_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  dd <- d
  dd$total_cn[dd$total_cn == 0] <- 150000
  dd$total_cn[is.na(dd$total_cn)] <- 0
  dd <- as(dd, "GRanges")
  cov <- coverage(dd, weight = dd$total_cn)
  dd1 <- as(cov, "GRanges")
  dd1 <- as.data.frame(dd1)
  dd1 <- dd1[dd1$score != 0, ]
  dd1 <- dd1[, c(1, 2, 3, 6)]
  names(dd1) <- names(d)[1:4]
  dd1$total_cn[dd1$total_cn == 150000] <- 0
  cnv <- dd1
  
  CN_data <- CNVsegs(chrom = as.character(cnv$chromosome),
                     start = cnv$start,
                     end   = cnv$end,
                     total_cn = cnv$total_cn)
  
  # Run ShatterSeek (genome version hg38)
  chromothripsis <- shatterseek(SV.sample = SV_data, seg.sample = CN_data, genome = "hg38")
  
  # Write output
  write.table(chromothripsis@chromSummary,
              file = paste0(out_prefix, ".out.txt"),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  saveRDS(chromothripsis, file = paste0(out_prefix, ".out.rds"))
}


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
sample_list_file <- args[1]
output_basename <- ifelse(length(args) >= 2, args[2], "test_shatterseek")

# Read sample information
sample_info <- read.table(sample_list_file, sep = "\t", header = FALSE,
                          stringsAsFactors = FALSE, comment.char = "")
colnames(sample_info) <- c("sample_id", "sv_file", "cnv_file")

output_dir <- "result"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  svf  <- sample_info$sv_file[i]
  cnvf <- sample_info$cnv_file[i]
  out_pref <- file.path(output_dir, sid)
  cat(sprintf("  Processing sample %s ...\n", sid))
  run_shatterseek_sample(svf, cnvf, out_pref)
}


all_files <- list.files(output_dir, pattern = "\\.out\\.txt$", full.names = TRUE)
if (length(all_files) == 0) {
  stop("No .out.txt files found. Please check that ShatterSeek ran successfully.")
}

# Initialize merged dataframe using header from first file
first_file <- all_files[1]
merged_df <- read.table(first_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
first_sample <- gsub("\\.out\\.txt$", "", basename(first_file))
merged_df <- cbind(sample = first_sample, merged_df)

# Append remaining files, skip headers, filter rows with insufficient fields (NF>10 condition)
for (f in all_files[-1]) {
  df <- read.table(f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  sample_name <- gsub("\\.out\\.txt$", "", basename(f))
  df <- df[apply(df, 1, function(x) sum(!is.na(x)) > 10), ]
  if (nrow(df) > 0) {
    df <- cbind(sample = sample_name, df)
    merged_df <- rbind(merged_df, df)
  }
}

# Write merged results
merged_file <- paste0(output_basename, "_result.txt")
write.table(merged_df, file = merged_file, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)

shattercall <- read.table(merged_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract p-value columns
ft   <- shattercall$pval_fragment_joins
ent  <- shattercall$chr_breakpoint_enrichment
ext  <- shattercall$pval_exp_chr
exct <- shattercall$pval_exp_cluster

# FDR adjustment
shattercall$ft_fdr   <- p.adjust(ft,   method = "fdr")
shattercall$ent_fdr  <- p.adjust(ent,  method = "fdr")
shattercall$ext_fdr  <- p.adjust(ext,  method = "fdr")
shattercall$exct_fdr <- p.adjust(exct, method = "fdr")

# Calculate total intra-chromosomal rearrangements
shattercall$total_intra <- shattercall$number_DEL + shattercall$number_DUP + shattercall$number_h2hINV + shattercall$number_t2tINV

# Write final result
final_file <- paste0(output_basename, "_result_fdr.txt")
write.table(shattercall, file = final_file, sep = "\t", quote = FALSE, row.names = FALSE)

