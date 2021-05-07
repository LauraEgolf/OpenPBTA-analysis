### Set root & analysis directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")
plots_dir <- file.path(analysis_dir, "plots", "05-candidate-regions/")

### Read in ShatterSeek results & list of chromoth objects
chromoth_combined <- read.table(file.path(analysis_dir, "results", "shatterseek_results_per_chromosome.txt"), 
                                  head=T, sep="\t", stringsAsFactors = F) 

chromoth_obj_list <- readRDS(file = file.path(root_dir, "scratch", "chromoth_obj_list.rds"))

### Define custom plotting function
plotChromothripsisCall <- function(confidence, chromoth_calls, i){
  bioid_to_plot <- chromoth_calls[i, "Kids_First_Biospecimen_ID"]
  chr_to_plot <- chromoth_calls[i, "chrom"]
  chromothripsis <- chromoth_obj_list[[bioid_to_plot]]
  plots <- plot_chromothripsis(ShatterSeek_output = chromothripsis, chr = chr_to_plot)
  plot_combined <- gridExtra::arrangeGrob( plots[[1]],
                                           plots[[2]],
                                           plots[[3]],
                                           plots[[4]], nrow=4,ncol=1,heights=c(0.2,.4,.4,.4))
  ggsave(paste0(plots_dir, confidence, "_", bioid_to_plot, "_chr", chr_to_plot, ".pdf"), plot_combined)
}

### Plot the first 10 high-confidence calls
chromoth_calls_hc <- subset(chromoth_combined, call_high_conf==1)
for (i in 1:10){
  plotChromothripsisCall("high_conf", chromoth_calls_hc, i)
}

### Plot the first 10 low-confidence calls
chromoth_calls_lc <- subset(chromoth_combined, call_low_conf==1)
for (i in 1:10){
  plotChromothripsisCall("low_conf", chromoth_calls_lc, i)
}

### Note that the ideograms are actually *hg19*
