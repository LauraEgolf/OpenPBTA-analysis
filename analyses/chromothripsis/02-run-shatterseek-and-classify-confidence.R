# Laura Egolf (@LauraEgolf) 2021
# Partially adapted from Yang Yang (@yangyangclover) 2020 
# https://github.com/AlexsLemonade/OpenPBTA-analysis/blob/master/analyses/sv-analysis/02-shatterseek.R

### Install ShatterSeek (use code from Dockerfile)
# BiocManager::install("graph")
# R -e "withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS='true'), remotes::install_github('parklab/ShatterSeek', ref = '83ab3effaf9589cc391ecc2ac45a6eaf578b5046', dependencies = TRUE))"


## ===================== Load Packages =====================
#library(devtools)
library(ShatterSeek)
library(tidyverse)   ### add in functions so I don't have to load entire library


## ===================== Root Directory =====================
# Detect the ".git" folder -- this will in the project root directory.
# Use this as the root directory to ensure proper sourcing of functions no
# matter where this is called from
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")

## ===================== Load  Independent Specimen List =====================
independent_specimen_list <- read.table(file.path(root_dir, "data", "independent-specimens.wgs.primary-plus.tsv"), 
                                        header = TRUE, sep = "\t", stringsAsFactors = F)
# bioid including all sample's names will be used later
bioid <- unique(independent_specimen_list$Kids_First_Biospecimen_ID)


## ===================== Load CNV File =====================
# Read cnv consensus file
cnvconsensus <- read_tsv(file.path(root_dir, "data", "pbta-cnv-consensus.seg.gz"))

# Shatterseek cannot work with NA copy number, remove rows with NA copy number
cnvconsensus <- cnvconsensus[!is.na(cnvconsensus$copy.num),]

# Choose independent specimens 
cnvconsensus <- cnvconsensus[cnvconsensus$ID %in% bioid,]

# Subset bioid to only samples that have CNV data
  # Note that 20 samples are not included in the CNV consensus because they failed QC 
  # for 2 callers (see analyses/copy_number_consensus_call/results/uncalled_samples.tsv).
  # So this analysis includes 777 samples instead of the full 797 in the independent
  # specimens list.
bioid <- bioid[bioid %in% cnvconsensus$ID]

# Remove chrY (ShatterSeek does not recognize)
cnvconsensus <- cnvconsensus[cnvconsensus$chrom != "chrY",]
# Remove "chr" notation (ShatterSeek does not recognize)
cnvconsensus$chrom <- gsub("chr","",cnvconsensus$chrom)

# # Check whether adjacent CNV segments have different CN values (requirement for ShatterSeek)
# checkCNV_logical <- logical()
# num <- dim(cnvconsensus)[1]-1
# for ( i in 1:num) {
#   checkCNV_logical[i] <-
#     (cnvconsensus[i, "loc.end"] == cnvconsensus[i+1, "loc.start"]) &
#     (cnvconsensus[i, "copy.num"] == cnvconsensus[i+1, "copy.num"]) &
#     (cnvconsensus[i, "chrom"] == cnvconsensus[i+1, "chrom"])
# }
# # TRUE for 16 segment pairs - should probably merge these...
# temp <- which(checkCNV_logical)
# temp2 <- c(temp, temp+1)
# temp3 <- unique(temp2[order(temp2)])
# View(cnvconsensus[temp3,])
# # They're all CN=2 regions - not sure if this matters


## ===================== Run ShatterSeek, combine and output results =====================
total <- length(bioid) # Total number of samples to run
count <- 0 # Keep track of sample count
chromoth_combined <- data.frame() # Merge results from different samples into one df
chromoth_obj_list <- vector("list", total) # Store shatterseek output as list of chromoth objects
names(chromoth_obj_list) <- bioid

# Loop through sample list and run ShatterSeek
for (i in bioid) {
  count=count+1
  print(paste0("Running: ", i, " (", count, " of ", total, ")"))
  
  # Read SV file for current sample
  sv_shatterseek <- read.table(file.path(root_dir, "scratch","sv-vcf",paste(i,"_withoutYandM.tsv",sep="")),sep="\t",header=TRUE)
  
  #  Subset CNV dataframe to current sample
  cnv_shatterseek <-  cnvconsensus[cnvconsensus$ID == i,]
  
  # If CNV or SV file is empty, jump into next loop
  if (nrow(cnv_shatterseek) == 0 | nrow(sv_shatterseek) == 0) {
      print(paste0(i," is missing CNV or SV data"))
    next;
  }

  # Build SV and CNV objects
  SV_data <-
    SVs(
      chrom1 = as.character(sv_shatterseek$chrom1),
      pos1 = as.numeric(sv_shatterseek$pos1),
      chrom2 = as.character(sv_shatterseek$chrom2),
      pos2 = as.numeric(sv_shatterseek$pos2),
      SVtype = as.character(sv_shatterseek$SVtype),
      strand1 = as.character(sv_shatterseek$strand1),
      strand2 = as.character(sv_shatterseek$strand2)
    )
  CN_data <-
    CNVsegs(
      chrom = as.character(cnv_shatterseek$chrom),
      start = cnv_shatterseek$loc.start,
      end = cnv_shatterseek$loc.end,
      total_cn = cnv_shatterseek$copy.num
    )
  
  # Run shatterseek
  chromothripsis <- shatterseek(SV.sample=SV_data,seg.sample=CN_data)

  # Save summary in a combined dataframe
  chromoth_summary <- chromothripsis@chromSummary
  chromoth_summary$Kids_First_Biospecimen_ID <- i
  chromoth_combined <- rbind(chromoth_combined, chromoth_summary)
  
  # Save full chromoth object as a named list
  chromoth_obj_list[[count]] <- chromothripsis
  
}

# Note: ShatterSeek only reports one chromothripsis region per chromosome. Also, the output (chromoth_combined)
# contains one row per chromosome per sample, even though not all chromosomes have a chromothripsis region. 

# Remove newline characters from column 'inter_other_chroms_coords_all'
# (I think they were inserted for plotting purposes, but they make it hard to read & write the data)
chromoth_combined$inter_other_chroms_coords_all <- gsub("\n", ";", chromoth_combined$inter_other_chroms_coords_all)

# Save list of chromoth objects in scratch directory
saveRDS(chromoth_obj_list, file = file.path(root_dir, "scratch", "chromoth_obj_list.rds"))


## ===================== Classify high and low confidence chromothripsis regions =====================
# Cutoffs used here are described briefly in ShatterSeek tutorial and in more detail 
# in Cortes-Ciriano et al (Supplemental Note).
# There is one set of cutoffs for low confidence chromothripsis regions, and two different 
# sets of criteria for high confidence regions.

### Add FDR correction
chromoth_combined$fdr_fragment_joins <- p.adjust(chromoth_combined$pval_fragment_joins, method = "fdr")
chromoth_combined$fdr_chr_breakpoint_enrichment <- p.adjust(chromoth_combined$chr_breakpoint_enrichment, method = "fdr")
chromoth_combined$fdr_exp_cluster <- p.adjust(chromoth_combined$pval_exp_cluster, method = "fdr")

### Define logic vector for each cutoff

## Low Confidence Cutoff
LC_cutoff <- 
  # At least 6 interleaved intrachromosomal SVs:
  ((chromoth_combined$clusterSize_including_TRA - chromoth_combined$number_TRA) >= 6) &
      # Note: this is equivalent to adding DEL, DUP, h2hINV, and t2tINV
  # At least 4 adjacent segments oscillating between 2 CN states
  (chromoth_combined$max_number_oscillating_CN_segments_2_states >= 4) & 
  # Not significant for the fragment joins test (even distribution of SV types)
  (chromoth_combined$fdr_fragment_joins > 0.2) & 
  # Significant for either the chromosomal enrichment or the exponential distribution of breakpoints test:
  (chromoth_combined$fdr_chr_breakpoint_enrichment < 0.2 | chromoth_combined$fdr_exp_cluster < 0.2)
      # Note: Use pval_exp_cluster, not pval_exp_chr

## High Confidence Cutoff 1
## Note: Same as low confidence cutoff, but more stringent for oscillating CN states
HC_cutoff1 <- 
  # At least 6 interleaved intrachromosomal SVs:
  ((chromoth_combined$clusterSize_including_TRA - chromoth_combined$number_TRA) >= 6) &
  # At least 7 adjacent segments oscillating between 2 CN states:
  (chromoth_combined$max_number_oscillating_CN_segments_2_states >= 7) & 
  # Not significant for the fragment joins test (even distribution of SV types)
  (chromoth_combined$fdr_fragment_joins > 0.2) & 
  # Significant for either the chromosomal enrichment or the exponential distribution of breakpoints test:
  (chromoth_combined$fdr_chr_breakpoint_enrichment < 0.2 | chromoth_combined$fdr_exp_cluster < 0.2)

## High Confidence Cutoff 2
HC_cutoff2 <-
  # At least 3 interleaved intrachromosomal SVs and at least 4 interchromosomal SVs:
  ((chromoth_combined$clusterSize_including_TRA - chromoth_combined$number_TRA) >= 3 & 
     chromoth_combined$number_TRA >= 4) &
  # At least 7 adjacent segments oscillating between 2 CN states:
  (chromoth_combined$max_number_oscillating_CN_segments_2_states >= 7) & 
  # Not significant for the fragment joins test (even distribution of SV types)
  (chromoth_combined$fdr_fragment_joins > 0.2)  

### Annotate each row of ShatterSeek results dataframe with chromothripsis call based on cutoffs for 
### high confidence, low confidence cutoff, or all confidence 
# Note "low_conf" reports calls that surpass low-confidence threshold but *not* high-confidence threshold
chromoth_combined$call_all_conf <- 0
chromoth_combined[which(LC_cutoff | HC_cutoff1 | HC_cutoff2), "call_all_conf"] <- 1
chromoth_combined$call_high_conf <- 0
chromoth_combined[which(HC_cutoff1 | HC_cutoff2), "call_high_conf"] <- 1
chromoth_combined$call_low_conf <- 0
chromoth_combined[which(LC_cutoff & !(HC_cutoff1 | HC_cutoff2)), "call_low_conf"] <- 1

### Create new dataframe with sample-level summary
# For each confidence call set: 
# Count the number of chromothripsis regions per sample
# Create a logical variable indicating whether or not each sample has >=1 chromothripsis region

chromoth_per_sample <- chromoth_combined %>% 
  dplyr::group_by(Kids_First_Biospecimen_ID) %>%
  dplyr::summarize_at(c("call_all_conf", "call_high_conf", "call_low_conf"), sum)

names(chromoth_per_sample) <- c("Kids_First_Biospecimen_ID", 
                                "count_regions_all_conf", "count_regions_high_conf", "count_regions_low_conf")

chromoth_per_sample$any_regions_all_conf <- chromoth_per_sample$count_regions_all_conf>0
chromoth_per_sample$any_regions_high_conf <- chromoth_per_sample$count_regions_high_conf>0
chromoth_per_sample$any_regions_low_conf <- chromoth_per_sample$count_regions_low_conf>0

# Merge "any_regions_*" into one column containing info on each confidence level
chromoth_per_sample$any_regions_merged <- "No Calls"
chromoth_per_sample[chromoth_per_sample$any_regions_low_conf, "any_regions_merged"] <- "Low Confidence"
chromoth_per_sample[chromoth_per_sample$any_regions_high_conf, "any_regions_merged"] <- "High Confidence"


## =====================  Write out results =====================

# ShatterSeek results with chromothripsis call defined above
write.table(chromoth_combined, file.path(analysis_dir, "results", "shatterseek_results_per_chromosome.txt"), 
            sep="\t", quote=F, row.names=F)

# Per-sample summary of chromothripsis calls
write.table(chromoth_per_sample, file.path(analysis_dir, "results", "chromothripsis_summary_per_sample.txt"), 
            sep="\t", quote=F, row.names=F)
