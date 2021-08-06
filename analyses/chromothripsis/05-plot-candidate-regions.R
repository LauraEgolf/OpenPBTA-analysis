### Load packages
library(ShatterSeek)
library(gridExtra)

### Define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")
plots_dir <- file.path(analysis_dir, "plots", "05-candidate-regions/")

### Read in chromothripsis regions and subset to low- or high-confidence calls
chromoth_combined <- read.table(file.path(analysis_dir, "results", "shatterseek_results_per_chromosome.txt"), 
                                head=T, sep="\t", stringsAsFactors = F) 
chromoth_calls <- subset(chromoth_combined, call_any_conf==1)

### Read in list of Shatterseek chromoth objects - need these for ShatterSeek::plot_chromothripsis()
chromoth_obj_list <- readRDS(file = file.path(root_dir, "scratch", "chromoth_obj_list.rds"))

### Define custom function - just a wrapper for ShatterSeek::plot_chromothripsis()
plotChromothripsisCall <- function(bioid_to_plot, chr_to_plot, confidence){
  if(is.na(confidence)){
    print("Warning: Call does not meet confidence threshold")
  }
  chromothripsis <- chromoth_obj_list[[bioid_to_plot]]
  plots <- plot_chromothripsis(ShatterSeek_output = chromothripsis, chr = chr_to_plot)
  plot_combined <- arrangeGrob(plots[[1]],
                               plots[[2]],
                               plots[[3]],
                               plots[[4]], nrow=4, ncol=1, heights=c(0.2,.4,.4,.4))
  ggplot2::ggsave(paste0(plots_dir, bioid_to_plot, "_chr", chr_to_plot, "_", confidence, ".pdf"), plot_combined)
}

### Select 10 random samples that have chromothripsis calls, and plot the candidate region(s) for each
set.seed(1234)
bioid_examples <- sample(unique(chromoth_calls$Kids_First_Biospecimen_ID), 10)
chromoth_calls_examples <- subset(chromoth_calls, Kids_First_Biospecimen_ID %in% bioid_examples)
number_of_examples <- dim(chromoth_calls_examples)[1]
for (i in 1:number_of_examples){
  bioid_to_plot <- chromoth_calls_examples[i, "Kids_First_Biospecimen_ID"]
  chr_to_plot <- chromoth_calls_examples[i, "chrom"]
  confidence <- ifelse(chromoth_calls_examples[i, "call_high_conf"] == 1, "high_conf", 
                 ifelse(chromoth_calls_examples[i, "call_low_conf"] == 1, "low_conf", 
                 NA))
  plotChromothripsisCall(bioid_to_plot, chr_to_plot, confidence)
}
