### Set root & analysis directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")
plots_dir <- file.path(analysis_dir, "plots", "breakpoint-data")

### Define Magrittr pipe
`%>%` <- dplyr::`%>%`

### Read in metadata_chromoth
# TO DO 
# Decide whether to write this out or repeat code


### Read in CNV and SV breakpoint data

breakpoint_data_dir <- file.path(root_dir, "analyses", "chromosomal-instability", "breakpoint-data")
cnv_densities <- read.table(file.path(breakpoint_data_dir, "cnv_breaks_densities.tsv"),
                            head=T, sep="\t", stringsAsFactors = F)
sv_densities <- read.table(file.path(breakpoint_data_dir, "sv_breaks_densities.tsv"),
                           head=T, sep="\t", stringsAsFactors = F)


### Merge breakpoint data into metadata/chromothripsis dataframe

# Reformat column names before merging
names(cnv_densities)[1] <- "Kids_First_Biospecimen_ID"
names(cnv_densities)[4] <- "cnv_breaks_count"
names(sv_densities)[1] <- "Kids_First_Biospecimen_ID"
names(sv_densities)[4] <- "sv_breaks_count"

merge <- metadata_chromoth %>% 
  dplyr::inner_join(cnv_densities, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::inner_join(sv_densities, by = "Kids_First_Biospecimen_ID") 


### Plot the number of chromothripsis regions per sample along with CNV breaks and/or SV breaks

# Scatterplot - chromothripsis regions, CNV breaks, SV breaks
p <- merge %>%
  dplyr::arrange(count_regions_all_conf) %>%
  # Sort by number of chromothripsis calls so that the samples with chromothripsis are more visible
  ggplot2::ggplot(ggplot2::aes(x = cnv_breaks_count, y = sv_breaks_count, color = as.factor(count_regions_all_conf))) +
  ggplot2::geom_point(shape=1, size=2) + 
  ggplot2::ylim(c(0, 5000)) +
  ggplot2::scale_color_brewer(palette="YlOrRd", name="# Chromothripsis\nRegions") +
  ggplot2::theme_light() +
  ggplot2::xlab("# CNV Breaks") + 
  ggplot2::ylab("# SV Breaks")
ggplot2::ggsave(file.path(plots_dir, "count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf"), p)

# Stripchart - chromothripsis regions, CNV breaks
p <- merge %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(count_regions_all_conf), 
                               y = cnv_breaks_count, color = as.factor(count_regions_all_conf))) +
  ggplot2::geom_jitter() +
  ggplot2::scale_color_brewer(palette="YlOrRd") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab("# Chromothripsis Regions") + 
  ggplot2::ylab("# CNV Breaks") 
ggplot2::ggsave(file.path(plots_dir, "count_chromothripsis_cnv_breaks_stripchart.pdf"), p)

# Stripchart - chromothripsis regions, SV breaks
p <- merge %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(count_regions_all_conf), 
                               y = sv_breaks_count, color = as.factor(count_regions_all_conf))) +
  ggplot2::geom_jitter() +
  ggplot2::scale_color_brewer(palette="YlOrRd") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab("# Chromothripsis Regions") + 
  ggplot2::ylab("# SV Breaks") 
ggplot2::ggsave(file.path(plots_dir, "count_chromothripsis_sv_breaks_stripchart.pdf"), p)
