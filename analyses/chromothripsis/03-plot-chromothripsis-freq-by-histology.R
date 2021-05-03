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

### Plot proportion of samples with chromothripsis out of the total for each display group 
# Use low confidence events 

# Plot by display_group
p <- metadata_chromoth %>%
  dplyr::count(any_low_conf, display_group, hex_codes) %>%
  tidyr::pivot_wider(names_from = any_low_conf, values_from = n, values_fill=0) %>%
  dplyr::group_by(display_group, hex_codes) %>%
  dplyr::mutate(group_size = sum(`TRUE`, `FALSE`)) %>%
  dplyr::mutate(prop = `TRUE` / group_size) %>%
  dplyr::mutate(labels = paste0("n=", group_size)) %>%
  ggplot2::ggplot(ggplot2::aes(x = display_group, y = prop, fill = hex_codes)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_text(aes(label=labels), vjust=-0.2, size=2.5) + 
  ggplot2::scale_fill_identity() +
  ggplot2::ylab("Proportion of Tumors with Chromothripsis") +
  ggplot2::theme(axis.text.x = element_text(angle = 90))
ggsave("plots/chromothripsis_proportion_per_display_group_low_confidence.pdf", p)

# Use pathology_diagnosis instead (to compare to yangyangclover plots)
p <- metadata_chromoth %>%
  dplyr::count(any_low_conf, pathology_diagnosis, hex_codes) %>%
  tidyr::pivot_wider(names_from = any_low_conf, values_from = n, values_fill=0) %>%
  dplyr::group_by(pathology_diagnosis, hex_codes) %>%
  dplyr::mutate(group_size = sum(`TRUE`, `FALSE`)) %>%
  dplyr::filter(group_size >= 5) %>%   # Remove groups with <5
  dplyr::mutate(prop = `TRUE` / group_size) %>%
  dplyr::mutate(labels = paste0("n=", group_size)) %>%
  ggplot2::ggplot(ggplot2::aes(x = pathology_diagnosis, y = prop, fill = hex_codes)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::geom_text(aes(label=labels), vjust=-0.2, size=2.5) + 
  ggplot2::scale_fill_identity() +
  ggplot2::ylab("Proportion of Tumors with Chromothripsis") +
  ggplot2::theme(axis.text.x = element_text(angle = 90))
ggsave("plots/chromothripsis_proportion_per_pathology_diagnosis_low_confidence.pdf", p)

