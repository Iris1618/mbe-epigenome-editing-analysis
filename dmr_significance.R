# Cappelluti et al., Nature (2024), doi:10.1038/s41586-024-07087-8
# Read the three GEO coverage archives and the two source workbooks directly.
# In R: source("dmr_significance.R"); out <- run_dmr_significance("/path/to/data")
# Use source_dir if the XLSX files are in a different folder from the archives.

# Identify the GEO samples and the source table columns for each plot.
# Individual CpG coordinates are recovered from the data at run time.
regions <- data.frame(
  system = c(rep("ZFP-ETRs", 2), rep("EvoETR-8", 4)),
  region = paste0("DMR-", c(11, 14, 1, 2, 4, 5)),
  series = paste0("GSE", c(240935, 240935, 240933, 240934, 240934, 240933)),
  first_gsm = c(7712391, 7712391, 7712379, 7712385, 7712385, 7712379),
  amplicon = c("offTargets", "offTargets", "chr2", "chr18_tchh", "chr18_tchh", "chr18"),
  chromosome = c("chr11", "chr17", "chr2", "chr3", "chr18", "chr18"),
  sheet = c(rep("Ext Data Fig. 7", 2), rep("Fig. 5G", 4)),
  cells = c("H6:I44", "F6:G62", "C7:D14", "F7:G119", "L7:M16", "O7:P53"),
  allow_vehicle_mismatch = c(rep(FALSE, 5), TRUE)
)

read_region <- function(region, data_dir, extract_dir) {
  archive <- file.path(data_dir, paste0(region$series, "_RAW.tar"))
  members <- utils::untar(archive, list = TRUE, tar = "internal")
  sample_ids <- paste0("GSM", region$first_gsm + 0:5)
  files <- vapply(sample_ids, function(id) {
    suffix <- paste0("_", region$amplicon, ".bismark.cov.gz")
    matches <- members[startsWith(members, paste0(id, "_")) & endsWith(members, suffix)]
    if (length(matches) != 1L || basename(matches) != matches) {
      stop("Expected one file for ", id, " / ", region$amplicon, " in ", archive)
    }
    matches
  }, character(1))
  status <- utils::untar(archive, files = files, exdir = extract_dir, tar = "internal")
  if (status != 0L) stop("Could not extract ", archive)

  samples <- lapply(files, function(file) {
    connection <- gzfile(file.path(extract_dir, file), open = "rt")
    on.exit(close(connection))
    coverage <- utils::read.delim(connection, header = FALSE,
      col.names = c("chrom", "start", "end", "percent", "methylated", "unmethylated"))
    coverage <- coverage[coverage$chrom == region$chromosome, ]
    depth <- coverage$methylated + coverage$unmethylated
    if (anyNA(depth) || anyDuplicated(coverage$start)) stop("Invalid coverage in ", file)
    sites <- depth >= 10
    stats::setNames(100 * coverage$methylated[sites] / depth[sites], coverage$start[sites])
  })
  positions <- sort(as.integer(Reduce(intersect, lapply(samples, names))))
  percentages <- vapply(samples, function(sample) sample[as.character(positions)],
                        numeric(length(positions)))
  colnames(percentages) <- sample_ids
  list(positions = positions, percentages = percentages)
}

match_plotted_sites <- function(region, candidates, workbook) {
  published <- as.matrix(readxl::read_excel(workbook, sheet = region$sheet, range = region$cells,
                                           col_names = FALSE, .name_repair = "minimal"))
  if (!is.numeric(published) || anyNA(published)) {
    stop("Invalid plot values in ", region$sheet, " ", region$cells)
  }
  percentages <- candidates$percentages
  if (!length(candidates$positions)) stop("No shared covered sites for ", region$region)
  treatment <- rowMeans(percentages[, 1:3])
  vehicle <- rowMeans(percentages[, 4:6])
  tolerance_pp <- 0.002
  treatment_close <- abs(outer(published[, 1], treatment, `-`)) <= tolerance_pp
  vehicle_close <- abs(outer(published[, 2], vehicle, `-`)) <= tolerance_pp

  # Require one match in both groups. DMR-5 has a known vehicle discrepancy,
  # so allow a unique treatment match when the vehicle values do not agree.
  plot_rows <- seq_len(nrow(published))
  matched <- vapply(plot_rows, function(row) {
    hits <- which(treatment_close[row, ] & vehicle_close[row, ])
    if (!length(hits) && region$allow_vehicle_mismatch) {
      hits <- which(treatment_close[row, ])
    }
    if (length(hits) != 1L) stop("No unique match for ", region$region, " plot row ", row)
    hits
  }, integer(1))
  if (is.unsorted(matched, strictly = TRUE)) {
    stop("GEO positions do not match the source-table order for ", region$region)
  }
  audit <- data.frame(
    system = region$system, region = region$region, chromosome = region$chromosome,
    position = candidates$positions[matched], plot_row = plot_rows,
    treatment_geo_pct = treatment[matched], treatment_published_pct = published[, 1],
    vehicle_geo_pct = vehicle[matched], vehicle_published_pct = published[, 2],
    vehicle_matches = vehicle_close[cbind(plot_rows, matched)], row.names = NULL
  )
  list(percentages = percentages[matched, , drop = FALSE], audit = audit)
}

compare_region <- function(percentages, system, region) {
  # Average CpGs within each mouse before testing between independent mice.
  mouse_means <- colMeans(percentages)
  treatment <- mouse_means[1:3]
  vehicle <- mouse_means[4:6]
  difference <- mean(treatment) - mean(vehicle)
  welch <- stats::t.test(treatment, vehicle, alternative = "two.sided", var.equal = FALSE)

  # Sensitivity analysis assuming labels are exchangeable under the null.
  permutations <- utils::combn(6L, 3L, function(i) {
    mean(mouse_means[i]) - mean(mouse_means[-i])
  })
  permutation_p <- mean(abs(permutations) >= abs(difference) - 1e-10)

  list(
    summary = data.frame(system, region, n_cpg = nrow(percentages),
      treatment_pct = mean(treatment), vehicle_pct = mean(vehicle),
      difference_pp = difference, ci95_low_pp = welch$conf.int[1],
      ci95_high_pp = welch$conf.int[2], p_welch = welch$p.value,
      p_exact_permutation = permutation_p),
    per_mouse = data.frame(system, region,
      group = rep(c("Treatment", "Vehicle"), each = 3L),
      sample = names(mouse_means), methylation_pct = unname(mouse_means))
  )
}

find_source_workbook <- function(source_dir, filenames) {
  files <- file.path(source_dir, filenames)
  found <- files[file.exists(files)]
  if (!length(found)) stop("Missing source workbook ", filenames[1], " in ", source_dir)
  found[1]
}

run_dmr_significance <- function(data_dir, source_dir = data_dir) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop('Install the Excel reader first with install.packages("readxl").')
  }
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  archives <- paste0(unique(regions$series), "_RAW.tar")
  missing <- archives[!file.exists(file.path(data_dir, archives))]
  if (length(missing)) stop("Missing GEO archive(s): ", paste(missing, collapse = ", "))

  fig5 <- find_source_workbook(source_dir, c("41586_2024_7087_MOESM15_ESM.xlsx", "SourceDataFig5.xlsx"))
  ed7 <- find_source_workbook(source_dir, c("41586_2024_7087_MOESM22_ESM.xlsx", "SourceDataExtendedFig7.xlsx"))
  pcsk9 <- readxl::read_excel(fig5, sheet = "Fig. 5E", range = "B5:K49",
                            col_names = FALSE, .name_repair = "minimal")
  if (anyNA(pcsk9) || anyDuplicated(pcsk9[[1]]) ||
      !identical(range(pcsk9[[1]]), c(106463706, 106464363))) {
    stop("Unexpected CpG data in Fig. 5E")
  }
  pcsk9 <- as.matrix(pcsk9[, 2:10])
  colnames(pcsk9) <- paste0(rep(c("ZFP", "Evo", "Vehicle"), each = 3L), "_col", 1:3)

  temp_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  extract_dir <- tempfile("dmr_coverage_", tmpdir = temp_root)
  if (!dir.create(extract_dir)) stop("Could not create a temporary extraction directory")
  extract_dir <- normalizePath(extract_dir, winslash = "/", mustWork = TRUE)
  if (dirname(extract_dir) != temp_root) stop("Unexpected temporary extraction directory")
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)
  off_target <- lapply(seq_len(nrow(regions)), function(i) {
    region <- regions[i, ]
    coverage <- read_region(region, data_dir, extract_dir)
    workbook <- if (region$sheet == "Fig. 5G") fig5 else ed7
    match_plotted_sites(region, coverage, workbook)
  })
  off_tests <- lapply(seq_len(nrow(regions)), function(i) {
    compare_region(off_target[[i]]$percentages, regions$system[i], regions$region[i])
  })
  zfp_pcsk9 <- compare_region(pcsk9[, c(1:3, 7:9)], "ZFP-ETRs", "DMR-3/4/5 (Pcsk9)")
  evo_pcsk9 <- compare_region(pcsk9[, c(4:6, 7:9)], "EvoETR-8", "DMR-3 (Pcsk9)")
  tests <- c(list(zfp_pcsk9), off_tests[1:4], list(evo_pcsk9), off_tests[5:6])
  results <- do.call(rbind, lapply(tests, `[[`, "summary"))
  results$p_holm <- stats::p.adjust(results$p_welch, method = "holm")
  results$significant_holm <- results$p_holm < 0.05
  per_mouse <- do.call(rbind, lapply(tests, `[[`, "per_mouse"))
  site_matches <- do.call(rbind, lapply(off_target, `[[`, "audit"))

  # Repeat DMR-5 without the three positions with discordant vehicle means.
  dmr5_index <- which(regions$allow_vehicle_mismatch)
  dmr5 <- off_target[[dmr5_index]]
  keep <- dmr5$audit$vehicle_matches
  sensitivity <- compare_region(dmr5$percentages[keep, , drop = FALSE],
                                regions$system[dmr5_index], regions$region[dmr5_index])$summary
  sensitivity$p_holm <- tail(stats::p.adjust(c(head(results$p_welch, -1), sensitivity$p_welch),
                                           method = "holm"), 1)
  sensitivity$significant_holm <- sensitivity$p_holm < 0.05

  print(results, row.names = FALSE, digits = 5)
  cat("\nDMR-5 sensitivity analysis (", sum(keep), " concordant CpGs)\n", sep = "")
  print(sensitivity, row.names = FALSE, digits = 5)
  invisible(list(results = results, per_mouse = per_mouse, site_matches = site_matches,
                 dmr5_sensitivity = sensitivity))
}
out <- run_dmr_significance(
  data_dir = ".",
  source_dir = "C:/Users/Iriss/Documents/ChatGPT/MBE_OralPre/tmp/offtarget_reanalysis"
)
View(out$results)