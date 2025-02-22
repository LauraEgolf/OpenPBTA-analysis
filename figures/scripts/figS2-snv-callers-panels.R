# S. Spielman for CCDL 2022
#
# Makes pdf panels for supplementary Figure S2, specifically those that are derived from the `snv-callers` analysis module.

library(tidyverse)


# Directories -------------------------------------------------------------------
# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
output_dir <- file.path(root_dir, "figures", "pdfs", "supp", "figs2", "panels")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Data directory
data_dir <- file.path(root_dir, "data")

# Analysis directores
analyses_dir <- file.path(root_dir, "analyses")
snv_callers_dir <- file.path(analyses_dir, "snv-callers")

# Scratch directory
scratch_dir <- file.path(root_dir, "scratch", "snv-callers")

# Palette directory
palette_dir <- file.path(root_dir, "figures", "palettes")


## Define final PDF files ------------------------------------------------------
pbta_upset_pdf <- file.path(output_dir, "pbta_upset_plot.pdf")
tcga_upset_pdf <- file.path(output_dir, "tcga_upset_plot.pdf")
pbta_vaf_cor_matrix_pdf <- file.path(output_dir, "pbta_vaf_cor_matrix.pdf")
tcga_vaf_cor_matrix_pdf <- file.path(output_dir, "tcga_vaf_cor_matrix.pdf")
pbta_vaf_distribution_plot_pdf <- file.path(output_dir, "pbta_vaf_distribution_plot.pdf")
tcga_vaf_distribution_plot_pdf <- file.path(output_dir, "tcga_vaf_distribution_plot.pdf")
lancet_wxs_wgs_plot_pdf <- file.path(output_dir, "lancet_wxs_wgs_plot.pdf")


## Read in data and load items --------------------------------------------------------------

# Read in clinical data
metadata <- read_tsv(file.path(data_dir, "pbta-histologies.tsv"),
                     guess_max = 10000)


# Binary palette
binary_palette <- read_tsv(file.path(palette_dir,
                                     "binary_color_palette.tsv"))


# Define columns to join by during data processing
join_cols <- c("Chromosome",
              "Start_Position",
              "Reference_Allele",
              "Allele",
              "Tumor_Sample_Barcode",
              "Variant_Classification")

# Loop over the two datasets
for (dataset in c("tcga", "pbta")) {
  print(paste0("===================",dataset,"===================="))
  # Set up DB connection ------------------------------------------------
  db_file <- "snv_db.sqlite"
  if (dataset == "tcga") {
    db_file <- paste0("tcga_", db_file)
  }
  con <- DBI::dbConnect(RSQLite::SQLite(),
                        file.path(scratch_dir, db_file))


  # Load and process database tables ------------------------------------
  print("Loading database tables")
  # Read in database tables, retaining only columns we need
  strelka <- tbl(con, "strelka") %>%
    select(join_cols, "VAF")

  lancet <- tbl(con, "lancet") %>%
    select(join_cols, "VAF")

  mutect <- tbl(con, "mutect") %>%
    select(join_cols, "VAF")

  if (dataset == "pbta") {

    # only PBTA has vardict
    vardict <- tbl(con, "vardict") %>%
      select(join_cols, "VAF")

    # Source a script that will full join these and create `all_caller`
    print("Joining database tables")
    source(file.path(snv_callers_dir, "util", "full_join_callers.R"))

    # data frame and add column to keep track of the original rows/mutations
    all_caller_df <- all_caller %>%
      as.data.frame() %>%
      rowid_to_column("index")

  } else if (dataset == "tcga") {
    # Need to `as.data.frame()` the TCGA ones since they are small enough to work with directly
    all_caller_df <- as.data.frame(strelka) %>%
      dplyr::full_join(as.data.frame(mutect), by = join_cols,
                       suffix = c("_strelka", "_mutect")) %>%
      dplyr::full_join(as.data.frame(lancet),
                       by = join_cols) %>%
      dplyr::rename(VAF_lancet = VAF) %>%
      # We'll use this to keep track of the original rows/mutations
      tibble::rowid_to_column("index")

  }


  ## Upset plots -------------------------------------------------------
  print("Preparing data for upset plot")
  detect_mat <- all_caller_df %>%
    # Bring over VAF columns and the index
    select(starts_with("VAF_")) %>%
    as.matrix()

  # Store the indices as dimnames
  dimnames(detect_mat)[[1]] <- all_caller_df$index

  # Turn into detected or not
  detect_mat <- !is.na(detect_mat)


  # Plot from `detect_mat`
  print("Making upset plot")
  # Set up a list how UpSetR wants it
  detect_list <- list(
    lancet = which(detect_mat[, "VAF_lancet"]),
    mutect = which(detect_mat[, "VAF_mutect"]),
    strelka = which(detect_mat[, "VAF_strelka"])
  )

  # add vardict if pbta, and define a shared plot_file
  if (dataset == "pbta") {
    # include vardict
    detect_list[["vardict"]] <- which(detect_mat[, "VAF_vardict"])
    plot_file <- pbta_upset_pdf
  } else if (dataset == "tcga") {
    plot_file <- tcga_upset_pdf
  }

  # Make plot and export to file
  # don't use the print() statement from original notebook - this adds a blank page in the DPF
  pdf(plot_file, width = 10, height = 5)
  UpSetR::upset(
    UpSetR::fromList(detect_list),
    order.by = "freq",
    text.scale = 1.2,
    point.size = 4,
    mainbar.y.label = ""
  )
  dev.off()



  ## VAF distribution plots ------------------------------------------
  print("Making VAF distribution plot")
  vaf_df <- all_caller_df %>%
    select(index, starts_with("VAF_"), Variant_Classification) %>%
    # Make long format
    gather(key = "caller", value = "vaf", -index, -Variant_Classification) %>%
    # Drop the `VAF_`
    mutate(caller = gsub("VAF_", "", caller)) %>%
    drop_na(vaf)

  vaf_plot <- ggplot(vaf_df) +
    aes(x = caller, y = vaf) +
    geom_violin(fill = "gray70") + # just a fill
    xlab("SNV Caller") +
    ylab("VAF") +
    ggpubr::theme_pubr()

  # export
  if (dataset == "pbta") {
    ggsave(pbta_vaf_distribution_plot_pdf, vaf_plot, width = 6, height = 4)
  } else if (dataset == "tcga") {
    ggsave(tcga_vaf_distribution_plot_pdf, vaf_plot, width = 6, height = 4)
  }


  ## VAF correlation plots --------------------------------------------
  print("Making VAF correlation plot")
  # Correlate VAFs across callers pbta
  cor_vaf <- all_caller_df %>%
    select(starts_with("VAF_")) %>%
    # 1:n() column for spread/gather bookkeeping which we do to remove `VAF_` from strips in this plot
    mutate(n = 1:n()) %>%
    # we want to remove
    gather(caller, vaf, -n) %>%
    mutate(caller = str_replace_all(caller, "VAF_", "")) %>%
    spread(caller, vaf) %>%
    select(-n) %>%
    GGally::ggpairs(aes(alpha = 0.05)) +
    ggpubr::theme_pubr() +
    cowplot::panel_border() +
    # slightly smaller text to avoid label overlap
    theme(axis.text = element_text(size = rel(0.8)))

  # export
  if (dataset == "pbta") {
    # try saving as TIFF with 300 DPI
	  print("saving TIFF")
    ggsave(file.path(output_dir, "pbta_vaf_cor_matrix.tiff"),
	   cor_vaf,
	   width = 9, height = 6, dpi = 300,
	   compression = "lzw")
    #ggsave(pbta_vaf_cor_matrix_pdf, cor_vaf, width = 6, height = 4)
  } else if (dataset == "tcga") {
    ggsave(tcga_vaf_cor_matrix_pdf, cor_vaf, width = 9, height = 6)
  }


  ## Lancet WXS/WGS plot for PBTA data only ------------------------------------------
  if (dataset == "pbta") {
    print("Making lancet WXS/WGS plot")
    # Retrieve all the participant IDs for participants that have both WGS and WXS data.
    matched_participants <- metadata %>%
      filter(experimental_strategy != "RNA-Seq") %>%
      group_by(Kids_First_Participant_ID) %>%
      summarize(strategies = paste0(experimental_strategy, collapse = ",")) %>%
      filter(grepl("WXS", strategies) & grepl("WGS", strategies)) %>%
      pull(Kids_First_Participant_ID)

    # Get the biospecimen IDS for these participants.
    biospecimens <- metadata %>%
      filter(Kids_First_Participant_ID %in% matched_participants) %>%
      pull(Kids_First_Biospecimen_ID)

    # Set up the Lancet data from the SQL database and only keep the biospecimens we identified.
    lancet <- tbl(con, "lancet") %>%
      select(
        join_cols, "VAF" #matches `cols_to_keep` in original notebook
      ) %>%
      inner_join(
        select(
          tbl(con, "samples"),
          Tumor_Sample_Barcode = Kids_First_Biospecimen_ID,
          experimental_strategy,
          short_histology,
          Kids_First_Participant_ID
          )
      ) %>%
      filter(Tumor_Sample_Barcode %in% biospecimens) %>%
      as.data.frame()

    # Participant and sample counts
    n_participants <- length(unique(lancet$Kids_First_Participant_ID))
    n_samples <- length(unique(lancet$Tumor_Sample_Barcode))

    # Prep binary colors
    colors <- binary_palette$hex_codes[binary_palette$color_names != "na_color"]

    # Plot and export
    vaf_plot <- ggplot(lancet) +
      aes(
        x = experimental_strategy,
        y = VAF,
        fill = experimental_strategy
      ) +
      geom_violin() +
      scale_fill_manual(values = colors) +
      labs(
        title = "Lancet participants with WGS and WXS",
        subtitle = glue::glue("{n_samples} samples from {n_participants} patients"),
        x = "Experimental Strategy",
        y = "VAF") +
      ggpubr::theme_pubr() +
      theme(legend.position = "none")

    ggsave(lancet_wxs_wgs_plot_pdf,
           vaf_plot,
           width = 5, height = 6)
  }


}









