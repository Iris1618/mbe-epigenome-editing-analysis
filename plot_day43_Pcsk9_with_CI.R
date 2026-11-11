# Locate this script when sourced or run from another working directory.
script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_files <- Filter(function(x) is.character(x) && length(x) == 1L,
                       lapply(sys.frames(), function(frame) frame$ofile))
script_file <- if (length(source_files)) tail(source_files, 1L)[[1L]] else
  if (length(script_args)) sub("^--file=", "", script_args[1L]) else NA_character_
script_dir <- if (!is.na(script_file))
  dirname(normalizePath(script_file, mustWork = TRUE)) else getwd()
data_dirs <- unique(c(script_dir, getwd(), file.path(getwd(), "analysis")))
data_dirs <- data_dirs[file.exists(file.path(data_dirs, "SourceDataFig5.xlsx"))]
if (!length(data_dirs)) stop("Cannot find SourceDataFig5.xlsx beside the script or in analysis/.")
analysis_dir <- normalizePath(data_dirs[1L], mustWork = TRUE)
local_library <- file.path(analysis_dir, "r-library")
if (dir.exists(local_library)) .libPaths(c(local_library, .libPaths()))

library(readxl)
library(ggplot2)

# 1. Day 43: C10:H10 = ZFP-ETRs; I10:N10 = EvoETR-8.
# Already vehicle-normalized; multiply by 100 to get percentages.
raw <- read_excel(file.path(analysis_dir, "SourceDataFig5.xlsx"), sheet = "Fig. 5B",
                  range = "C10:N10", col_names = FALSE, .name_repair = "minimal")
zfp <- as.numeric(unlist(raw[1, 1:6])) * 100
evo <- as.numeric(unlist(raw[1, 7:12])) * 100
groups <- c("ZFP-ETRs", "EvoETR-8")
dat <- data.frame(group = factor(rep(groups, each = 6), levels = groups),
                  pcsk9 = c(zfp, evo))
means <- data.frame(group = factor(groups, levels = groups),
                    mean = c(mean(zfp), mean(evo)), sd = c(sd(zfp), sd(evo)))

# 2. Same exploratory unpaired, two-sided Welch comparison.
# Positive difference = additional suppression with EvoETR-8.
# The CI does not propagate uncertainty in the shared vehicle denominator.
test <- t.test(zfp, evo, paired = FALSE, var.equal = FALSE)
effect <- mean(zfp) - mean(evo)
lo <- test$conf.int[1]
hi <- test$conf.int[2]

effect_label <- sprintf("%.1f pp (95%% CI %.1f to %.1f)", effect, lo, hi)
cat(effect_label, "\n")
cat(sprintf("Exploratory two-sided Welch test: p = %.4f (unadjusted)\n", test$p.value))

# Keep the original significance star and show the rounded p value.
p_label <- if (test$p.value < 0.05)
  sprintf("* p = %.3f < 0.05", test$p.value) else
  sprintf("p = %.3f", test$p.value)

# 3. Shared colours, fonts and export settings.
# Headings/captions are placed separately as editable text on the slide.
navy <- "#062E52"
purple <- "#80499B"
blue <- "#3A68A0"
font_family <- "Arial"
export_dpi <- 450
theme_slide <- theme_classic(base_size = 17, base_family = font_family) +
  theme(axis.title = element_text(colour = navy, size = 17),
        axis.text = element_text(colour = navy, size = 16),
        axis.line = element_line(colour = navy, linewidth = 0.55),
        axis.ticks = element_line(colour = navy, linewidth = 0.5),
        legend.position = "none",
        plot.margin = margin(8, 12, 8, 8),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))

# 4. Individual mice, with mean +/- SD.
# Keep the original deterministic horizontal jitter; no vertical jitter.
p_individual <- ggplot(dat, aes(group, pcsk9, colour = group)) +
  geom_errorbar(data = means,
                aes(x = group, ymin = mean - sd, ymax = mean + sd),
                inherit.aes = FALSE, width = 0.10, linewidth = 0.75, colour = navy) +
  geom_point(position = position_jitter(width = 0.11, height = 0, seed = 2026),
             size = 4, alpha = 1) +
  geom_errorbar(data = means, aes(x = group, ymin = mean, ymax = mean),
                inherit.aes = FALSE, width = 0.45, linewidth = 1.1, colour = navy) +
  annotate("segment", x = 1, xend = 2, y = 86, yend = 86,
           colour = navy, linewidth = 0.65) +
  annotate("segment", x = c(1, 2), xend = c(1, 2), y = 83, yend = 86,
           colour = navy, linewidth = 0.65) +
  annotate("text", x = 1.5, y = 94, label = p_label,
           colour = navy, size = 6, fontface = "bold", family = font_family) +
  scale_colour_manual(values = c("ZFP-ETRs" = blue, "EvoETR-8" = purple)) +
  scale_x_discrete(limits = groups,
                   labels = c("ZFP-ETRs" = "ZFP-ETRs\n(n = 6)",
                              "EvoETR-8" = "EvoETR-8\n(n = 6)"),
                   expand = expansion(add = 0.55)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Remaining plasma PCSK9\n(% of vehicle)") +
  theme_slide + theme(axis.title.y = element_text(margin = margin(r = 12)))

# 5. Moderately compact horizontal estimate and 95% confidence interval.
# Keep the estimate and CI visible without restoring the tall empty panel.
effect_annotation <- sprintf("%.1f percentage points\n95%% CI: %.1f to %.1f",
                             effect, lo, hi)
p_effect <- ggplot(data.frame(estimate = effect, lower = lo, upper = hi), aes(y = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.65, colour = "#A5AFB9") +
  geom_segment(aes(x = lower, xend = upper, yend = 1), colour = purple, linewidth = 1.25) +
  geom_segment(aes(x = lower, xend = lower, y = 0.975, yend = 1.025),
               colour = purple, linewidth = 0.9) +
  geom_segment(aes(x = upper, xend = upper, y = 0.975, yend = 1.025),
               colour = purple, linewidth = 0.9) +
  geom_point(aes(x = estimate), colour = purple, size = 4.5) +
  annotate("text", x = effect, y = 1.23, label = effect_annotation,
           colour = navy, size = 4.8, fontface = "bold", family = font_family,
           lineheight = 1.1) +
  scale_x_continuous(limits = c(-10, 50), breaks = seq(-10, 50, 10),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = c(0.87, 1.34), breaks = NULL, expand = c(0, 0)) +
  labs(x = "Additional suppression (percentage points)", y = NULL) +
  theme_slide +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 17, margin = margin(t = 10)),
        plot.margin = margin(6, 18, 6, 18))

# 6. Display in RStudio and export directly from ggplot2.
# Keep the original filenames; give the CI figure a moderate height.
if (interactive()) {
  print(p_individual)
  print(p_effect)
}
ggsave(file.path(analysis_dir, "day43_individual_mice.png"), p_individual,
       width = 5.4, height = 4.7, dpi = export_dpi, bg = "white", type = "cairo")
ggsave(file.path(analysis_dir, "day43_difference_95CI.png"), p_effect,
       width = 6.2, height = 2.6, dpi = export_dpi, bg = "white", type = "cairo")
ggsave(file.path(analysis_dir, "day43_individual_mice.pdf"), p_individual,
       width = 5.4, height = 4.7, bg = "white", device = cairo_pdf)
ggsave(file.path(analysis_dir, "day43_difference_95CI.pdf"), p_effect,
       width = 6.2, height = 2.6, bg = "white", device = cairo_pdf)
cat("P-value label:", p_label, "\n")
cat("Figures saved in:", analysis_dir, "\n")
