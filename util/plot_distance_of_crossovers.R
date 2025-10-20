#!/usr/bin/env Rscript

# ============================================================
#  crossover_distance_histogram.R
#  Analyse distances between consecutive crossovers
#  within phase groups from HaplotypeTools output
#
#  Usage:
#    Rscript crossover_distance_histogram.R path/to/your_file.tab
#
#  Expected columns (tab-delimited):
#    1. phase_group   (arbitrary group ID)
#    2. contig        (chromosome / scaffold)
#    3. position      (integer base position)
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a file name, e.g. Rscript crossover_distance_histogram.R crossovers.tab")
}
file <- args[1]

# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tools)
})

# ---- read file ----
dat <- fread(file, select = 1:3, header = FALSE)
setnames(dat, c("phase_group", "contig", "position"))

# Make sure position is numeric
dat[, position := as.numeric(position)]

# Drop NAs or malformed rows
dat <- dat[!is.na(position) & position > 0]

# ---- compute distances ----
setorder(dat, phase_group, contig, position)

# compute distance to previous site within same phase group + contig
dat[, dist_to_prev := position - shift(position, type = "lag"), by = .(phase_group, contig)]

# keep only positive distances (ignore 0 or NA)
distances <- dat[!is.na(dist_to_prev) & dist_to_prev > 0, dist_to_prev]

# summary
cat("Total crossovers:", nrow(dat), "\n")
cat("Valid distances:", length(distances), "\n")
cat("Summary (bp):\n")
print(summary(distances))
cat("Quantiles (bp):\n")
print(quantile(distances, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)))

# ---- output file ----
pdf_file <- paste0(file_path_sans_ext(file), "_distances.pdf")

# ---- plot ----
p <- ggplot(data.frame(distance = distances), aes(x = distance)) +
  geom_histogram(bins = 100, color = "black", fill = "steelblue") +
  scale_x_continuous(trans = "log10",
                     breaks = 10^seq(0, ceiling(log10(max(distances, na.rm=TRUE))), by = 1),
                     labels = scales::comma) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of distances between consecutive crossovers",
    x = "Distance between crossovers (bp, log10 scale)",
    y = "Count"
  )

pdf(pdf_file, width = 8, height = 5)
print(p)
dev.off()

cat("Histogram saved to:", pdf_file, "\n")
