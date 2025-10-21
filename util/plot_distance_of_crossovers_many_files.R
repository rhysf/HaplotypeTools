#!/usr/bin/env Rscript

# ============================================================
#  crossover_distance_histogram_multi.R
#  Combine multiple HaplotypeTools crossover files
#  and plot the aggregated distance distribution
#
#  Usage:
#    Rscript crossover_distance_histogram_multi.R file_list.txt
#
#  where 'file_list.txt' is a plain text file with one
#  tab- or newline-separated path to each file, e.g.:
#    sample1.tab
#    sample2.tab
#    sample3.tab
#
#  Each input file should have:
#    1. phase_group
#    2. contig
#    3. position
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a list file, e.g. Rscript crossover_distance_histogram_multi.R list_of_files.txt")
}
list_file <- args[1]

# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(tools)
})

# read list of files
files <- fread(list_file, header = FALSE)$V1
if (length(files) == 0) stop("No files found in the list!")

# ---- helper: compute distances per file ----
compute_distances <- function(file) {
  dat <- tryCatch(fread(file, select = 1:3, header = FALSE), error = function(e) NULL)
  if (is.null(dat)) {
    warning(paste("Skipping unreadable file:", file))
    return(NULL)
  }
  setnames(dat, c("phase_group", "contig", "position"))
  dat[, position := as.numeric(position)]
  dat <- dat[!is.na(position) & position > 0]
  setorder(dat, phase_group, contig, position)
  dat[, dist_to_prev := position - shift(position, type = "lag"), by = .(phase_group, contig)]
  dat[!is.na(dist_to_prev) & dist_to_prev > 0, .(distance = dist_to_prev)]
}

# gather all distances
all_distances <- rbindlist(lapply(files, compute_distances), use.names = TRUE, fill = TRUE)

if (nrow(all_distances) == 0) stop("No valid distances found in any file.")

# ---- summary ----
cat("Number of files processed:", length(files), "\n")
cat("Total distances combined:", nrow(all_distances), "\n")
cat("Summary (bp):\n")
print(summary(all_distances$distance))
cat("Quantiles (bp):\n")
print(quantile(all_distances$distance, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)))

#output
pdf_file <- paste0(file_path_sans_ext(basename(list_file)), "_combined_distances.pdf")

# ---- plot ----
p <- ggplot(all_distances, aes(x = distance)) +
  geom_histogram(bins = 100, color = "black", fill = "steelblue") +
  scale_x_continuous(
    trans = "log10",
    breaks = 10^seq(0, ceiling(log10(max(all_distances$distance, na.rm = TRUE))), by = 1),
    labels = scales::comma
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("Distribution of crossover distances (", length(files), " samples combined)"),
    x = "Distance between consecutive crossovers (bp, log10 scale)",
    y = "Count"
  )

pdf(pdf_file, width = 8, height = 5)
print(p)
dev.off()

cat("Histogram saved to:", pdf_file, "\n")