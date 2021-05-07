### Set root & analysis directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "chromothripsis")
plots_dir <- file.path(analysis_dir, "plots", "06-circos")


# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Source custom circos plot functions
source(file.path(root_dir, "analyses", "chromosomal-instability", "util", "circos-plots.R"))

# Read in the metadata
metadata <- readr::read_tsv(file.path(data_dir, "pbta-histologies.tsv"))

# Read in the CNV data
cnv_df <- data.table::fread(
  file.path(data_dir, 
            "pbta-cnv-consensus.seg.gz"),
  data.table = FALSE
) 

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

### Read, filter, and reformat SV data
sv_df <- data.table::fread(
  file.path(data_dir, "pbta-sv-manta.tsv.gz"),
  data.table = FALSE
)

sv_df <- subset(sv_df, FILTER=="PASS")

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


# Chromothripsis regions meeting low confidence threshold
chromoth_calls <- subset(chromoth_combined, call_low_conf==1)

samples_for_examples <- chromoth_calls[1:5, "Kids_First_Biospecimen_ID"]

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

cnv_df <- cnv_df %>%
  # Make column that specifies the color for each value
  dplyr::mutate(num_color_key = color_fun(copy.num))

# Let's determine how many levels this factor column has
n_levels <- length(levels(cnv_df$status))

# Set up a palette based on number of factor levels
palette_col <- RColorBrewer::brewer.pal(n_levels, name = "Accent")

# Let's make a key to recode by based on levels
palette_key <- palette_col

# Have the factor levels be the names
names(palette_key) <- levels(cnv_df$status)

cnv_df <- cnv_df %>%
  # Make column that specifies the color for each factor level
  dplyr::mutate(fac_color_key = dplyr::recode(status, !!!palette_key))




set.seed(2020)
samples_for_examples <- sample(cnv_df$ID, 5)


circos_map_plot(
  df = cnv_df,
  add_track = FALSE,
  samples_col = "ID",
  sample_names = samples_for_examples[1:3], # What samples we are plotting.
  chr_col = "chrom",
  start_col = "loc.start",
  end_col = "loc.end",
  y_val = "copy.num",
  track_height = .15,
  type = "rect", # Changed this to rect
  rect_height = .2, # Optionally can change height with this argument. Default is +_ 0.4
  color_col = "num_color_key", 
  cytoband = FALSE # Turning off the cytoband here. Default is TRUE
)

circos_map_transloc(transloc_df,
                    add_track = TRUE, # We change this to true to add on to our already existing plot
                    sample_names = samples_for_examples[3],
                    samples_col = "biospecimen_id1",
                    chr_col_1 = "chrom1", # Need to specify which column is the first and second location for each
                    chr_col_2 = "chrom2",
                    start_col_1 = "start1",
                    start_col_2 = "start2",
                    end_col_1 = "end1",
                    end_col_2 = "end2"
)
