# MBE EvoETR Reanalysis

This repository contains three exploratory R analyses of EvoETR-8 and ZFP-ETRs using published mouse data for our MBE oral presentation.

## Paper

Cappelluti, M. A. et al. [Durable and efficient gene silencing in vivo by hit-and-run epigenome editing](https://doi.org/10.1038/s41586-024-07087-8). *Nature* **627**, 416-423 (2024).

The study asked whether transient delivery of engineered transcriptional repressors could durably silence *Pcsk9* in mouse liver. The authors reported sustained suppression of circulating PCSK9 after a single dose and observed that suppression and promoter methylation persisted after liver regeneration. They also developed EvoETR-8, an all-in-one editor, and compared it with the original ZFP-ETR combination.

## Analyses

### Day 43 plasma PCSK9

[R script](plot_day43_Pcsk9_with_CI.R) · Figure 5b · Six mice per group

We extracted the day 43 plasma PCSK9 values from the published workbook. The values were already normalized to vehicle, so we multiplied them by 100 to express them as percentages. We used a two-sided, unpaired Welch t-test and subtracted the EvoETR-8 mean from the ZFP-ETR mean. A positive difference represents additional PCSK9 suppression with EvoETR-8.

**Result.** EvoETR-8 produced an estimated **21.4 percentage points** of additional suppression (95% CI 2.1-40.8, unadjusted p = 0.036). The interval treats the shared vehicle denominator as fixed.

### Day 1 plasma ALT

[R script](plot_day1_ALT_with_CI.R) · Extended Data Figure 5a · Six mice per group

We selected day 1 after inspecting the reported ALT time course. We compared the two groups using a two-sided, unpaired Welch t-test on log-transformed alanine aminotransferase (ALT). We exponentiated the difference and its confidence limits to obtain the EvoETR-8-to-ZFP-ETR geometric mean ratio and its 95% confidence interval.

**Result.** The estimated ratio was **2.28** (95% CI 0.58–8.97, unadjusted p = 0.207). The interval includes 1 and is wide, leaving the size of the day 1 difference uncertain.

Both plasma comparisons are exploratory and use one time point per outcome.

### In vivo methylation at selected DMRs

[R script](dmr_significance.R) · Figure 5g and Extended Data Figure 7 · Three mice per group

For each region, we averaged methylation across the plotted CpGs within each mouse. We then compared the treatment group with its vehicle control using a two-sided, unpaired Welch t-test. The analysis covers six off-target regions and two on-target comparisons at *Pcsk9*.

**Result.** At an unadjusted p < 0.05, EvoETR-8 showed higher methylation at three of four tested off-target regions (DMR-2, DMR-4, and DMR-5). Neither tested ZFP-ETR off-target region reached that threshold. The presentation uses these exploratory, unadjusted results. After Holm correction across all eight comparisons, two off-target regions still met the threshold (EvoETR-8 DMR-4 and DMR-5). The script also reports two-sided exact permutation p-values. The authors tested 2 of 18 ZFP-ETR candidates and 4 of 4 EvoETR-8 candidates in mice, which limits the comparison of in vivo specificity.

For EvoETR-8 DMR-5, three vehicle CpG means differed between the GEO coverage data and the published source table. The script uses the GEO values and repeats the analysis with these three sites excluded.

## Run the scripts

### Plasma PCSK9 and ALT

1. Install the R packages `readxl` and `ggplot2` with `install.packages(c("readxl", "ggplot2"))`. Figure export uses R's Cairo graphics support.
2. In the [original paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10937395/), go to **Source data** and download **Source Data Fig. 5** and **Source Data Extended Data Fig. 5**. Save them beside the scripts as `SourceDataFig5.xlsx` and `SourceDataExtendedFig5.xlsx`, respectively.
3. Run the scripts from the repository folder.

```sh
Rscript plot_day43_Pcsk9_with_CI.R
Rscript plot_day1_ALT_with_CI.R
```

Each script exports an individual-mouse plot and an effect estimate with a 95% confidence interval as PNG and PDF. The PCSK9 outputs are saved beside its workbook. The ALT script saves figures and CSV summaries in `figures/`.

### DMR methylation analysis

1. Install `readxl` if needed with `install.packages("readxl")` and place the following files in the repository folder. The GEO archives contain per-mouse Bismark methylation coverage files. Keep the archives intact.

   - [GSE240933_RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240933/suppl/GSE240933_RAW.tar) for EvoETR-8 DMR-1 and DMR-5.
   - [GSE240934_RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240934/suppl/GSE240934_RAW.tar) for EvoETR-8 DMR-2 and DMR-4.
   - [GSE240935_RAW.tar](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE240nnn/GSE240935/suppl/GSE240935_RAW.tar) for ZFP-ETRs DMR-11 and DMR-14.
   - **Source Data Fig. 5** and **Source Data Extended Data Fig. 7** from the [paper's Source data section](https://pmc.ncbi.nlm.nih.gov/articles/PMC10937395/). Keep the original workbook filenames or name them `SourceDataFig5.xlsx` and `SourceDataExtendedFig7.xlsx`, respectively.

2. Open `dmr_significance.R` in RStudio with the repository folder as the working directory. Set `source_dir = "."` in the call at the bottom of the script and click **Source**.

The results table reports treatment and vehicle means, differences in percentage points, 95% confidence intervals, and p-values. `out$per_mouse` contains the per-mouse values and `out$site_matches` records how the GEO CpGs matched the published plots.
