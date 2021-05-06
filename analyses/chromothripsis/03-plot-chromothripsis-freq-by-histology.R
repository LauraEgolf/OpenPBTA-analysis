### Set root & analysis directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")
plots_dir <- file.path(analysis_dir, "plots", "histology")

### Define Magrittr pipe
`%>%` <- dplyr::`%>%`

### Load packages
library(ggplot2)
library(RColorBrewer)

### Read in chromothripsis info per sample (# regions, chromothripsis yes/no)
chromoth_per_sample <- read.table(file.path(analysis_dir, "results", "chromothripsis_info_per_sample.txt"), 
                                  head=T, sep="\t", stringsAsFactors = F)

### Merge with standard color palettes and metadata

# Import standard color palettes & histology display groups for project
histology_label_mapping <- readr::read_tsv(
  file.path(root_dir, "figures", "palettes", "histology_label_color_table.tsv")) %>% 
  # Select just the columns we will need for plotting
  dplyr::select(Kids_First_Biospecimen_ID, display_group, display_order, hex_codes) %>% 
  # Reorder display_group based on display_order
  dplyr::mutate(display_group = forcats::fct_reorder(display_group, display_order))

# Merge with metadata file
metadata_file <- file.path(root_dir, "data", "pbta-histologies.tsv")
metadata <- readr::read_tsv(metadata_file, guess_max = 10000) %>%
  dplyr::inner_join(histology_label_mapping, by = "Kids_First_Biospecimen_ID")

# Append chromothripsis info to merged metadata file
metadata_chromoth <- dplyr::inner_join(metadata, chromoth_per_sample, by="Kids_First_Biospecimen_ID")


### Plot proportion of samples with chromothripsis out of the total for each histology 

### Plot by pathology_diagnosis (to compare to yangyangclover plots)
# Basic bar plot: all events (any confidence level)
p <- metadata_chromoth %>%
  dplyr::count(any_regions_all_conf, pathology_diagnosis, hex_codes) %>%
  tidyr::pivot_wider(names_from = any_regions_all_conf, values_from = n, values_fill=0) %>%
  dplyr::group_by(pathology_diagnosis, hex_codes) %>%
  dplyr::mutate(group_size = sum(`TRUE`, `FALSE`)) %>%
  dplyr::filter(group_size >= 5) %>%   # Remove groups with <5
  dplyr::mutate(prop = `TRUE` / group_size) %>%
  dplyr::mutate(labels = paste0(`TRUE`, " / ", group_size)) %>%
  ggplot(aes(x = pathology_diagnosis, y = prop, fill = hex_codes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=labels), vjust=-0.2, size=2.5) + 
  scale_fill_identity() +
  xlab(NULL) + 
  ylab("Proportion of Tumors with Chromothripsis") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95))
ggsave(file.path(plots_dir, "chromothripsis_proportion_per_pathology_diagnosis.pdf"), p)


### For the rest of the plots, plot by display_group


### Basic bar plot: all events (any confidence level)

p <- metadata_chromoth %>%
  dplyr::count(any_regions_all_conf, display_group, hex_codes) %>%
  tidyr::pivot_wider(names_from = any_regions_all_conf, values_from = n, values_fill=0) %>%
  dplyr::group_by(display_group, hex_codes) %>%
  dplyr::mutate(group_size = sum(`TRUE`, `FALSE`)) %>%
  dplyr::mutate(prop = `TRUE` / group_size) %>%
  dplyr::mutate(labels = paste0(`TRUE`, " / ", group_size)) %>%
  ggplot(aes(x = display_group, y = prop, fill = hex_codes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=labels), vjust=-0.2, size=2.5) + 
  scale_fill_identity() +
  xlab(NULL) + 
  ylab("Proportion of Tumors with Chromothripsis") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95))
ggsave(file.path(plots_dir, "chromothripsis_proportion_per_display_group.pdf"), p)


### Bar plot stacked by low vs. high confidence 

p <- ggplot(metadata_chromoth, aes(x = display_group, fill = hex_codes, 
                                   alpha = factor(any_regions_merged, 
                                                  levels=c("No Calls", "Low Confidence", "High Confidence")))) +
  geom_bar(position = "fill") +
  scale_fill_identity() +
  scale_alpha_manual(values=c(0, 0.5, 1), name="Confidence Level") +
  ylim(c(0,0.4)) +
  xlab(NULL) + 
  ylab("Proportion of Tumors") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95))
ggsave(file.path(plots_dir, "chromothripsis_proportion_per_display_group_withConfidence.pdf"), p)


### Bar plot stacked by number of chromothripsis regions

# Define color scale for # chromothripsis regions, but set "0" as transparent so the bar doesn't show
  # This may need to be updated if the ShatterSeek results change
colors <- brewer.pal(7, "YlOrRd")
colors[1] <- "#1C00ff00"

p <- ggplot(metadata_chromoth, aes(x = display_group, fill = as.factor(count_regions_all_conf))) +
  geom_bar(position = "fill") +
  scale_fill_manual(values=colors, name="# Chromothripsis\nRegions") +
  ylim(c(0,0.4)) +
  xlab(NULL) + 
  ylab("Proportion of Tumors") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, hjust=0.95))
ggsave(file.path(plots_dir, "chromothripsis_proportion_per_display_group_colorByCount.pdf"), p)
