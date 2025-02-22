## Chromothripsis Analysis (in progress)

**Module authors :**
Laura Egolf ([@LauraEgolf](https://github.com/LauraEgolf/)), Yang Yang ([@yangyangclover](https://github.com/yangyangclover))

This module runs ShatterSeek, identifies chromothripsis regions, and visualizes the results.

#### Inputs from data download
- `independent-specimens.wgs.primary-plus.tsv`
- `pbta-sv-manta.tsv`
- `pbta-cnv-consensus.seg`

#### Order of scripts in analysis
`01-process-sv-file.R` : This script reformats SV files for ShatterSeek input.

`02-run-shatterseek-and-classify-confidence.R` : This script runs ShatterSeek and classifies chromothripsis regions as high or low confidence based on criteria recommended by the authors. 

`03-plot-chromothripsis-by-histology.Rmd`: This script plots the proportion of tumors with chromothripsis across histology groups.

`04-plot-chromothripsis-and-breakpoint-data.Rmd`: This script plots the chromothripsis results alongside breakpoint data from `analyses/chromosomal-instability/`.

More scripts will be added for visualization.

#### Output tables
`results/shatterseek_results_per_chromosome.txt` : Each row reports ShatterSeek output for one chromosome from one sample. (Note that ShatterSeek only reports one candidate chromothripsis region per chromosome, and it reports a row for each chromosome even if that chromosome doesn't have a chromothripsis region.) Each row is also annotated with three columns indicating how this chromothripsis region was classified in `02-run-shatterseek-and-classify-confidence.R`:
- `call_any_conf` : Chromothripsis call of any confidence (0=no, 1=yes)
- `call_high_conf` : High-confidence chromothripsis call (0=no, 1=yes)
- `call_low_conf` : Low-confidence chromothripsis call (0=no, 1=yes) - these calls surpass low-confidence threshold but *not* high-confidence threshold

`results/chromothripsis_summary_per_sample.txt` : Each row reports the number of chromothripsis calls for one sample:
- `count_regions_any_conf` : Number of chromothripsis calls of any confidence
- `count_regions_high_conf` : Number of high-confidence chromothripsis calls
- `count_regions_low_conf` : Number of low-confidence chromothripsis calls
- `any_regions` : Summary variable indicating whether a sample has no calls, >=1 low-confidence call, or >=1 high-confidence call

#### Plots
`plots/03-histology/*`: Proportion of tumors with chromothripsis across histology groups, with a breakdown by confidence level (high or low confidence chromothripsis) and by the number of chromothripsis regions observed. (See also: `03-plot-chromothripsis-by-histology.html`)

`plots/04-breakpoint-data/*`: Scatterplots/stripcharts showing the number of chromothripsis regions per sample alongside the breakpoint counts calculated from SV and CNV data in `analyses/chromosomal-instability/`. (See also: `04-plot-chromothripsis-and-breakpoint-data.html`)
