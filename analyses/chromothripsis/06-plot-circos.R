### Adapted from @canssavy
# https://cansavvy.github.io/openpbta-notebook-concept/chromosomal-instability/01b-visualization-cnv-sv.nb.html

# Define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")
plots_dir <- file.path(analysis_dir, "plots", "06-circos/")

# Define Magrittr pipe
`%>%` <- dplyr::`%>%`

# Source custom circos plot functions
source(file.path(root_dir, "analyses", "chromosomal-instability", "util", "circos-plots.R"))

# Read in chromothripsis regions and subset to low- or high-confidence calls
chromoth_combined <- read.table(file.path(analysis_dir, "results", "shatterseek_results_per_chromosome.txt"), 
                                head=T, sep="\t", stringsAsFactors = F) 
chromoth_calls <- subset(chromoth_combined, call_any_conf==1)

# Read in the metadata
metadata <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"))

# Read in the CNV data
cnv_df <- data.table::fread(
  file.path(data_dir, 
            "pbta-cnv-consensus.seg.gz"),
  data.table = FALSE
) 

# Create a status variable based on copy number and tumor ploidy
cnv_df <- cnv_df %>% 
  # Only keep data that has values for copy.num
  dplyr::filter(!is.na(copy.num)) %>%
  # Tack on the ploidy and sex estimate information from metadata
  dplyr::inner_join(
    dplyr::select(metadata, 
                  Kids_First_Biospecimen_ID, 
                  tumor_ploidy,
                  germline_sex_estimate), 
    by = c("ID" = "Kids_First_Biospecimen_ID")) %>%
  # Create a status variable
  dplyr::mutate(status = dplyr::case_when(
    # when the copy number is less than inferred ploidy, mark this as a loss
    copy.num < tumor_ploidy ~ "loss",
    # if copy number is higher than ploidy, mark as a gain
    copy.num > tumor_ploidy ~ "gain",
    copy.num == tumor_ploidy ~ "neutral"
  ), 
  status = factor(status, levels = c("loss", "neutral", "gain")))

# Read and filter SV data
sv_df <- data.table::fread(
  file.path(data_dir, "pbta-sv-manta.tsv.gz"),
  data.table = FALSE
)
sv_df <- subset(sv_df, FILTER=="PASS")

# Make a translocation data.frame where both sets of coordinates for the translocation are in the same row.
transloc_df <- sv_df %>%
  dplyr::filter(SV.type == "BND") %>%
  dplyr::mutate(
    match_id = stringr::str_sub(ID, 0, -3),
    pair_num = stringr::str_sub(ID, -1)
  ) %>%
  dplyr::select(
    biospecimen_id = Kids.First.Biospecimen.ID.Tumor,
    chrom = SV.chrom,
    start = SV.start,
    end = SV.end,
    width = SV.length,
    match_id,
    pair_num
  )
transloc_df <- transloc_df %>%
  dplyr::filter(pair_num == 0) %>%
  dplyr::inner_join(dplyr::filter(transloc_df, pair_num == 1),
                    by = "match_id",
                    suffix = c("1", "2")
  )

# Make color palette based on 5 colors
palette_col <- RColorBrewer::brewer.pal(5, 
                                        name = "Accent" # Can change this palette 
                                        # Use RColorBrewer::display.brewer.all() to see options
)

# Make color ramp function based on quantiles of seg.mean and palette
color_fun <- circlize::colorRamp2(
  breaks = quantile(cnv_df$seg.mean, 
                    c(0.15, 0.35, 0.5, 0.65, 0.85), na.rm = TRUE), 
  colors = palette_col)

# Make column that specifies the color for each value
cnv_df <- cnv_df %>%
  dplyr::mutate(num_color_key = color_fun(copy.num))

# Let's determine how many levels this factor column has
n_levels <- length(levels(cnv_df$status))

# Set up a palette based on number of factor levels
palette_col <- RColorBrewer::brewer.pal(n_levels, name = "Accent")

# Let's make a key to recode by based on levels
palette_key <- palette_col

# Have the factor levels be the names
names(palette_key) <- levels(cnv_df$status)

# Make column that specifies the color for each factor level
cnv_df <- cnv_df %>%
  dplyr::mutate(fac_color_key = dplyr::recode(status, !!!palette_key))
# Convert to string (doesn't work as factor)
cnv_df$fac_color_key <- as.character(cnv_df$fac_color_key)

####### To Do: Shorten/clean up code above #######

# Make column that specifies color for high vs. low confidence chromothripsis calls (orange vs. yellow)
chromoth_calls$color_chromoth_regions <- NA
chromoth_calls[chromoth_calls$call_high_conf==1, "color_chromoth_regions"] <- "orange"
chromoth_calls[chromoth_calls$call_low_conf==1, "color_chromoth_regions"] <- "yellow"

# Make column with arbitrary y-value to plot chromothripsis call regions
chromoth_calls$y_val <- 1


### Define custom plotting function to build circos plot with chromothripsis regions,
### CNV segments, and translocation arcs

chromothripsisCircosPlot <- function(bioid_to_plot){
  # Initialize plot
  pdf(paste0(plots_dir, "circos_", bioid_to_plot, ".pdf"))
  
  # Plot chromothripsis regions
  circos_map_plot(
    df = chromoth_calls,
    add_track = FALSE,
    samples_col = "Kids_First_Biospecimen_ID",
    sample_names = bioid_to_plot, 
    chr_col = "chrom",
    start_col = "start",
    end_col = "end",
    y_val = "y_val",
    track_height = .05,
    type = "rect", 
    rect_height = .5, # Optionally can change height with this argument. Default is +_ 0.4
    color_col = "color_chromoth_regions", 
    cytoband = TRUE
  )
  
  # Add CNV segments
  circos_map_plot(
    df = cnv_df,
    add_track = TRUE,
    samples_col = "ID",
    sample_names = bioid_to_plot,
    chr_col = "chrom",
    start_col = "loc.start",
    end_col = "loc.end",
    y_val = "copy.num",
    track_height = .15,
    type = "rect", 
    rect_height = .2,
    color_col = "fac_color_key" 
  )
  
  # Add translocation arcs
  circos_map_transloc(transloc_df,
                      add_track = TRUE, 
                      sample_names = bioid_to_plot,
                      samples_col = "biospecimen_id1",
                      chr_col_1 = "chrom1", # Need to specify which column is the first and second location for each
                      chr_col_2 = "chrom2",
                      start_col_1 = "start1",
                      start_col_2 = "start2",
                      end_col_1 = "end1",
                      end_col_2 = "end2"
  )
  dev.off()
}

### Select 10 random samples that have chromothripsis calls and draw circos plot for each sample
set.seed(1234)
bioid_examples <- sample(unique(chromoth_calls$Kids_First_Biospecimen_ID), 10)
for (bioid in bioid_examples){
  chromothripsisCircosPlot(bioid)
}
