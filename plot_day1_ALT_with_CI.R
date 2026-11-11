# Day-1 ALT: individual mice and geometric-mean ratio with 95% CI.
# Source: Cappelluti et al., Nature 627, 416-423 (2024).
# https://doi.org/10.1038/s41586-024-07087-8
#
# Statistical scope: an exploratory, unpaired, two-sided Welch test on log(ALT).
# All 12 values are retained. No equivalence test or multiplicity adjustment.
# Day 1 was chosen after inspecting the reported time course.

script_directory <- function() {
  sourced_files <- lapply(sys.frames(), function(frame) frame$ofile)
  sourced_files <- Filter(function(x) is.character(x) && length(x) == 1L,
                          sourced_files)
  if (length(sourced_files)) {
    return(dirname(normalizePath(tail(sourced_files, 1L)[[1L]],
                                 winslash = "/", mustWork = TRUE)))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1L]),
                                 winslash = "/", mustWork = TRUE)))
  }
  getwd()
}

script_dir <- script_directory()
source_name <- "SourceDataExtendedFig5.xlsx"
data_candidates <- unique(c(
  file.path(script_dir, source_name),
  file.path(script_dir, "analysis", source_name),
  file.path(script_dir, "..", "..", "analysis", source_name),
  file.path(getwd(), source_name),
  file.path(getwd(), "analysis", source_name)
))
existing <- data_candidates[file.exists(data_candidates)]
if (!length(existing)) {
  stop("Cannot find ", source_name,
       ". Place it beside this script, then click Source again.")
}
source_file <- normalizePath(existing[1L], winslash = "/", mustWork = TRUE)

# Use the project's existing local package library when it is present.
local_library <- file.path(dirname(source_file), "r-library")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))
required <- c("readxl", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install the missing packages first: install.packages(c(",
       paste(sprintf('"%s"', missing), collapse = ", "), "))")
}
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(ggplot2))

# Register Arial for Windows/RStudio graphics as well as Cairo exports.
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("TT Arial"))
}

output_dir <- file.path(script_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Source cells: Ext Data Fig. 5A, D8:I8 = ZFP-ETRs; J8:O8 = EvoETR-8.
# These are day-1 plasma ALT measurements in U/L, six independent mice/group.
raw <- read_excel(source_file, sheet = "Ext Data Fig. 5A", range = "D8:O8",
                  col_names = FALSE, .name_repair = "minimal")
zfp <- as.numeric(unlist(raw[1, 1:6], use.names = FALSE))
evo <- as.numeric(unlist(raw[1, 7:12], use.names = FALSE))
stopifnot(length(zfp) == 6L, length(evo) == 6L,
          all(is.finite(c(zfp, evo))), all(c(zfp, evo) > 0))

# Exponentiating the Welch CI for mean log(ALT) difference produces the CI
# for the ratio of the two population geometric means. Ratio > 1 means
# higher ALT with EvoETR-8. It is not a ratio of toxicity or injury severity.
welch <- t.test(log(evo), log(zfp), paired = FALSE, var.equal = FALSE,
                alternative = "two.sided", conf.level = 0.95)
groups <- c("ZFP-ETRs", "EvoETR-8")
individual <- data.frame(group = factor(rep(groups, each = 6), levels = groups),
                         mouse_index = rep(1:6, 2), ALT_U_L = c(zfp, evo))
group_summary <- data.frame(group = factor(groups, levels = groups), n = c(6L, 6L),
                            geometric_mean_ALT_U_L = c(exp(mean(log(zfp))),
                                                       exp(mean(log(evo)))))
effect <- data.frame(
  ratio_evo_over_zfp = exp(mean(log(evo)) - mean(log(zfp))),
  ci95_lower = exp(welch$conf.int[1L]),
  ci95_upper = exp(welch$conf.int[2L]),
  p_two_sided_welch = welch$p.value,
  df_welch = unname(welch$parameter)
)

write.csv(individual, file.path(output_dir, "day1_ALT_individual_data.csv"),
          row.names = FALSE)
write.csv(group_summary, file.path(output_dir, "day1_ALT_group_summary.csv"),
          row.names = FALSE)
write.csv(effect, file.path(output_dir, "day1_ALT_effect.csv"), row.names = FALSE)

# Colours match the existing day-43 presentation slide.
navy <- "#062E52"
blue <- "#3A68A0"
purple <- "#80499B"
axis_text_size <- 30
axis_title_size <- 32
p_value_size <- 10.5

base <- theme_classic(base_size = 32, base_family = "Arial") +
  theme(axis.title = element_text(colour = navy, size = axis_title_size),
        axis.text = element_text(colour = navy, size = axis_text_size),
        axis.line = element_line(colour = navy, linewidth = 0.55),
        axis.ticks = element_line(colour = navy, linewidth = 0.5),
        legend.position = "none",
        plot.margin = margin(16, 24, 16, 12),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# Left panel: each point is one mouse; short horizontal lines are geometric
# means. Jitter affects horizontal position only and never changes ALT.
p_individual <- ggplot(individual, aes(group, ALT_U_L, colour = group)) +
  geom_point(position = position_jitter(width = 0.10, height = 0, seed = 2026),
             size = 4, alpha = 1) +
  geom_errorbar(data = group_summary,
                aes(x = group, ymin = geometric_mean_ALT_U_L,
                    ymax = geometric_mean_ALT_U_L),
                inherit.aes = FALSE, width = 0.43, linewidth = 1.1, colour = navy) +
  annotate("segment", x = 1, xend = 2, y = 25000, yend = 25000,
           colour = navy, linewidth = 0.7) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = 21000, yend = 25000,
           colour = navy, linewidth = 0.7) +
  annotate("text", x = 1.5, y = 42000,
           label = sprintf("p = %.3f %s 0.05", welch$p.value,
                           if (welch$p.value > 0.05) ">" else if (welch$p.value < 0.05) "<" else "="),
           colour = navy,
           size = p_value_size, fontface = "bold", family = "Arial") +
  scale_colour_manual(values = c("ZFP-ETRs" = blue, "EvoETR-8" = purple)) +
  scale_x_discrete(limits = groups,
                   labels = c("ZFP-ETRs\n(n = 6)", "EvoETR-8\n(n = 6)"),
                   expand = expansion(add = 0.35)) +
  scale_y_log10(limits = c(100, 80000), breaks = c(100, 1000, 10000),
                labels = c("100", "1,000", "10,000"),
                expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "ALT (U/L, log scale)") +
  base + theme(axis.title.y = element_text(margin = margin(r = 16)),
               axis.text.x = element_text(lineheight = 0.95,
                                          margin = margin(t = 10)))

# Right panel: log-ratio axis; 1 is the null value (equal geometric means).
# The interval crosses 1, so the difference is not statistically significant.
# Its width also shows that a substantial ALT increase cannot be excluded.
p_ci <- ggplot(effect) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.65,
             colour = "#A5AFB9") +
  geom_segment(aes(x = ci95_lower, xend = ci95_upper, y = 1, yend = 1),
               linewidth = 1.25, colour = purple) +
  geom_segment(aes(x = ci95_lower, xend = ci95_lower,
                   y = 0.95, yend = 1.05), linewidth = 0.9, colour = purple) +
  geom_segment(aes(x = ci95_upper, xend = ci95_upper,
                   y = 0.95, yend = 1.05), linewidth = 0.9, colour = purple) +
  geom_point(aes(x = ratio_evo_over_zfp, y = 1), size = 4.5, colour = purple) +
  scale_x_log10(limits = c(0.4, 10), breaks = c(0.5, 1, 2, 4, 8),
                labels = c("0.5", "1", "2", "4", "8"),
                expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0.86, 1.14), breaks = NULL,
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "Geometric mean ALT ratio\n(EvoETR-8 / ZFP-ETRs, log scale)", y = NULL) +
  base +
  theme(aspect.ratio = 0.20,
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = axis_title_size, lineheight = 0.95,
                                    margin = margin(t = 16)),
        axis.text.x = element_text(margin = margin(t = 8)),
        plot.margin = margin(16, 24, 16, 24))

# Slide-ready panels: no duplicate heading, allowing editable PPT headings.
ggsave(file.path(output_dir, "day1_ALT_individual.png"), p_individual,
       width = 6.8, height = 5.8, dpi = 450, bg = "white", type = "cairo")
ggsave(file.path(output_dir, "day1_ALT_ratio_ci.png"), p_ci,
       width = 12, height = 4.5, dpi = 450, bg = "white", type = "cairo")
ggsave(file.path(output_dir, "day1_ALT_individual.pdf"), p_individual,
       width = 6.8, height = 5.8, bg = "white", device = cairo_pdf)
ggsave(file.path(output_dir, "day1_ALT_ratio_ci.pdf"), p_ci,
       width = 12, height = 4.5, bg = "white", device = cairo_pdf)

# RStudio: two separate entries in the Plots pane history.
if (interactive()) {
  print(p_individual)
  print(p_ci)
}
