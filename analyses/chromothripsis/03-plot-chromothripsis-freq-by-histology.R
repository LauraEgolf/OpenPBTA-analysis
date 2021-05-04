# Need to add setup

### Read in chromothripsis info per sample (# regions, chromothripsis yes/no)
chromoth_per_sample <- read.table(file.path(analysis_dir, "results", "chromothripsis_info_per_sample.txt"), 
                                  head=T, sep="\t", stringsAsFactors = F)


### Setup for plotting

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
metadata <- subset(metadata, Kids_First_Biospecimen_ID %in% bioid)

# Append chromothripsis info to merged metadata file
metadata_chromoth <- dplyr::inner_join(metadata, chromoth_per_sample, by="Kids_First_Biospecimen_ID")

### Plot proportion of samples with chromothripsis (low confidence) out of the total for each display group 

# Plot by display_group
p <- metadata_chromoth %>%
  dplyr::count(any_low_conf, display_group, hex_codes) %>%
  tidyr::pivot_wider(names_from = any_low_conf, values_from = n, values_fill=0) %>%
  dplyr::group_by(display_group, hex_codes) %>%
  dplyr::mutate(group_size = sum(`TRUE`, `FALSE`)) %>%
  dplyr::mutate(prop = `TRUE` / group_size) %>%
  dplyr::mutate(labels = paste0(`TRUE`, " / ", group_size)) %>%
  ggplot2::ggplot(ggplot2::aes(x = display_group, y = prop, fill = hex_codes)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_text(aes(label=labels), vjust=-0.2, size=2.5) + 
  ggplot2::scale_fill_identity() +
  ggplot2::ylab("Proportion of Tumors with Chromothripsis") +
  ggplot2::theme_light() +
  ggplot2::theme(axis.text.x = element_text(angle = 90))
ggsave("plots/chromothripsis_proportion_per_display_group.pdf", p)

# Use pathology_diagnosis instead (to compare to yangyangclover plots)
p <- metadata_chromoth %>%
  dplyr::count(any_low_conf, pathology_diagnosis, hex_codes) %>%
  tidyr::pivot_wider(names_from = any_low_conf, values_from = n, values_fill=0) %>%
  dplyr::group_by(pathology_diagnosis, hex_codes) %>%
  dplyr::mutate(group_size = sum(`TRUE`, `FALSE`)) %>%
  dplyr::filter(group_size >= 5) %>%   # Remove groups with <5
  dplyr::mutate(prop = `TRUE` / group_size) %>%
  dplyr::mutate(labels = paste0(`TRUE`, " / ", group_size)) %>%
  ggplot2::ggplot(ggplot2::aes(x = pathology_diagnosis, y = prop, fill = hex_codes)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_text(aes(label=labels), vjust=-0.2, size=2.5) + 
  ggplot2::scale_fill_identity() +
  ggplot2::ylab("Proportion of Tumors with Chromothripsis") +
  ggplot2::theme_light() +
  ggplot2::theme(axis.text.x = element_text(angle = 90))
ggsave("plots/chromothripsis_proportion_per_pathology_diagnosis.pdf", p)



### Read in CNV and SV breakpoint data; reformat
breakpoint_data_dir <- file.path(root_dir, "analyses", "chromosomal-instability", "breakpoint-data")
cnv_densities <- read.table(file.path(breakpoint_data_dir, "cnv_breaks_densities.tsv"),
                                head=T, sep="\t", stringsAsFactors = F)
sv_densities <- read.table(file.path(breakpoint_data_dir, "sv_breaks_densities.tsv"),
                               head=T, sep="\t", stringsAsFactors = F)
names(cnv_densities)[1] <- "Kids_First_Biospecimen_ID"
names(cnv_densities)[4] <- "cnv_breaks_count"
names(sv_densities)[1] <- "Kids_First_Biospecimen_ID"
names(sv_densities)[4] <- "sv_breaks_count"

### Merge breakpoint data into metadata/chromothripsis dataframe
merge <- metadata_chromoth %>% 
  dplyr::inner_join(cnv_densities, by = "Kids_First_Biospecimen_ID") %>%
  dplyr::inner_join(sv_densities, by = "Kids_First_Biospecimen_ID") 

### Plot the number of chromothripsis regions (low confidence) per sample along with CNV breaks and/or SV breaks

# Scatterplot - chromothripsis regions, CNV breaks, SV breaks
p <- merge %>%
  dplyr::arrange(count_regions_low_conf) %>%
    # Sort by number of chromothripsis calls so that the samples with chromothripsis are more visible
  ggplot2::ggplot(ggplot2::aes(x = cnv_breaks_count, y = sv_breaks_count, color = as.factor(count_regions_low_conf))) +
  ggplot2::geom_point(shape=1, size=2) + 
  ggplot2::ylim(c(0, 5000)) +
  ggplot2::scale_color_brewer(palette="YlOrRd", name="Chromothripsis\nRegions Count") +
  ggplot2::theme_light() +
  ggplot2::xlab("CNV Breaks Count") + 
  ggplot2::ylab("SV Breaks Count")
ggsave("plots/count_chromothripsis_cnv_and_sv_breaks_scatterplot.pdf", p)

# Stripchart - chromothripsis regions, CNV breaks
p <- merge %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(count_regions_low_conf), 
                               y = cnv_breaks_count, color = as.factor(count_regions_low_conf))) +
  ggplot2::geom_jitter() +
  ggplot2::scale_color_brewer(palette="YlOrRd") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab("Chromothripsis Regions Count") + 
  ggplot2::ylab("CNV Breaks Count") 
ggsave("plots/count_chromothripsis_cnv_breaks_stripchart.pdf", p)

# Stripchart - chromothripsis regions, SV breaks
p <- merge %>%
  ggplot2::ggplot(ggplot2::aes(x = as.factor(count_regions_low_conf), 
                               y = sv_breaks_count, color = as.factor(count_regions_low_conf))) +
  ggplot2::geom_jitter() +
  ggplot2::scale_color_brewer(palette="YlOrRd") +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::xlab("Chromothripsis Regions Count") + 
  ggplot2::ylab("SV Breaks Count") 
ggsave("plots/count_chromothripsis_sv_breaks_stripchart.pdf", p)
  
