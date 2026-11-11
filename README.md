# MBE EvoETR Reanalysis

This repository contains two exploratory R analyses comparing EvoETR-8 and ZFP-ETRs using published mouse data for our MBE oral presentation.

## Paper

Cappelluti, M. A. et al. [Durable and efficient gene silencing in vivo by hit-and-run epigenome editing](https://doi.org/10.1038/s41586-024-07087-8). *Nature* **627**, 416–423 (2024).

The study asked whether transient delivery of engineered transcriptional repressors could durably silence *Pcsk9* in mouse liver. The authors reported sustained suppression of circulating PCSK9 after a single dose and observed that suppression and promoter methylation persisted after liver regeneration. They also developed EvoETR-8, an all-in-one editor, and compared it with the original ZFP-ETR combination.

## Analyses

### Day 43 plasma PCSK9

[R script](plot_day43_Pcsk9_with_CI.R) · Figure 5b · Six mice per group

We extracted the day 43 plasma PCSK9 values from the published workbook. The values were already normalized to vehicle, so we multiplied them by 100 to express them as percentages. We used a two-sided, unpaired Welch t-test and subtracted the EvoETR-8 mean from the ZFP-ETR mean. A positive difference represents additional PCSK9 suppression with EvoETR-8.

**Result.** EvoETR-8 produced an estimated **21.4 percentage points** of additional suppression (95% CI 2.1–40.8, unadjusted p = 0.036). The interval treats the shared vehicle denominator as fixed.

### Day 1 plasma ALT

[R script](plot_day1_ALT_with_CI.R) · Extended Data Figure 5a · Six mice per group

We selected day 1 after inspecting the reported ALT time course. We compared the two groups using a two-sided, unpaired Welch t-test on log-transformed alanine aminotransferase (ALT). We exponentiated the difference and its confidence limits to obtain the EvoETR-8-to-ZFP-ETR geometric mean ratio and its 95% confidence interval.

**Result.** The estimated ratio was **2.28** (95% CI 0.58–8.97, unadjusted p = 0.207). The interval includes 1 and is wide, leaving the size of the day 1 difference uncertain.

Both comparisons are exploratory and use one time point per outcome.

## Run the scripts

1. Install the R packages `readxl` and `ggplot2` with `install.packages(c("readxl", "ggplot2"))`. Figure export uses R's Cairo graphics support.
2. In the [original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10937395/), go to **Source data** and download **Source Data Fig. 5** and **Source Data Extended Data Fig. 5**. Save them beside the scripts as `SourceDataFig5.xlsx` and `SourceDataExtendedFig5.xlsx`, respectively.
3. Run the scripts from the repository folder.

```sh
Rscript plot_day43_Pcsk9_with_CI.R
Rscript plot_day1_ALT_with_CI.R
```

Each script exports an individual-mouse plot and an effect estimate with a 95% confidence interval as PNG and PDF. The PCSK9 outputs are saved beside its workbook. The ALT script saves figures and CSV summaries in `figures/`.
