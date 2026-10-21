# LUAD hypoxia-related nine-gene signature

R scripts for a retrospective study of a nine-gene score associated with survival in lung adenocarcinoma (LUAD). The analysis builds the score using TCGA-LUAD data, evaluates it in other public cohorts, and examines its biological context.

The scripts accompany:

**A Hypoxia-Related Nine-Gene Signature Associated With Overall Survival and Tumor Biology in Lung Adenocarcinoma: A Retrospective Observational Study**

The repository also includes single-cell, pathway, immune-context, mutation and predicted drug-response analyses, plus plotting support for the RT-qPCR hypoxia experiments. These are research analyses, not a clinical decision tool. The drug-response estimates are computational predictions, not measured patient treatment responses.

## Repository Contents

| Script | Purpose |
|---|---|
| `01_download_TCGA_data.R` | Download TCGA-LUAD expression and clinical data, and collect the Hallmark hypoxia gene list. |
| `02_differential_expression.R` | Compare tumor and normal expression with DESeq2 and plot the hypoxia-related genes. |
| `03_LASSO_Cox_model.R` | Screen genes for survival associations, fit the LASSO-Cox score, and plot TCGA survival. Includes the LDHA/GAPDH exclusion analysis. |
| `04_prognostic_evaluation.R` | Evaluate survival associations, calibration and discrimination (C-index); build a nomogram and clinical summary table. |
| `05_GSE68465_validation.R` | Apply TCGA-derived coefficients to GSE68465, using that cohort's own median to define risk groups. |
| `06_scRNA_analysis.R` | Examine the score across cell populations in LUAD single-cell data. |
| `07_pathway_enrichment.R` | Compare genes and pathways between risk groups using GO, KEGG, Hallmark GSEA and GSVA. |
| `08_immune_microenvironment.R` | Examine immune-related expression and ESTIMATE scores, and prepare/read CIBERSORTx and TIDE files. |
| `09_drug_sensitivity_qPCR.R` | Estimate drug responses with GDSC2/oncoPredict and plot RT-qPCR measurements. |
| `10_GEO_data_extraction.R` | Extract GSE68465 data from GEO and map expression probes to genes. |
| `11_model_benchmarking_TMB_GSE50081.R` | Run supporting model comparisons, TCGA mutation-burden/driver analyses and GSE50081 evaluation. |
| `12_make_supportive_figures.R` | Plot the supporting model, mutation and GSE50081 results. |
| `13_GSE50081_supportive_evaluation.R` | Run the separate supporting evaluation in GSE50081. |

Some filenames retain names from earlier development. Use the final manuscript and figure legends to identify which results belong to the submitted analysis.

## Data Sources

The main datasets are publicly available:

- TCGA-LUAD RNA-seq, clinical data, and masked somatic mutation data: Genomic Data Commons / TCGA.
- GSE68465, GSE50081, and GSE131907: Gene Expression Omnibus.
- Hallmark hypoxia gene set: MSigDB via `msigdbr`.
- GDSC2 expression and drug-response data: GDSC / oncoPredict resources.
- CIBERSORTx and TIDE outputs require running the corresponding external web tools using the input files exported by the scripts.

No new high-throughput sequencing data were generated in the study. RT-qPCR raw and processed data are provided in the manuscript supplementary materials.

## Dependencies

The scripts were developed in R. Main package dependencies include:

```text
TCGAbiolinks, SummarizedExperiment, DESeq2, limma, clusterProfiler,
enrichplot, msigdbr, org.Hs.eg.db, GSVA, survival, survminer, glmnet,
rms, timeROC, tableone, Seurat, oncoPredict, estimate,
GEOquery, hgu133a.db, ggplot2, ggpubr, ggrepel, ggvenn, pheatmap,
patchwork, dplyr, tidyr, tibble, stringr, scales, openxlsx, gridExtra
```

Package versions and changes to online databases can affect the output. Save `sessionInfo()` with your results so the R and package versions are recorded.

## Suggested Execution Order

1. Run `01_download_TCGA_data.R`.
2. Run `02_differential_expression.R`.
3. Run `03_LASSO_Cox_model.R`.
4. Run `04_prognostic_evaluation.R`.
5. Run `10_GEO_data_extraction.R`, then `05_GSE68465_validation.R`.
6. Run `06_scRNA_analysis.R` after obtaining the required GSE131907 matrix/annotation files.
7. Run `07_pathway_enrichment.R`.
8. Run `08_immune_microenvironment.R`; submit exported files to CIBERSORTx/TIDE where required, then rerun or continue with returned outputs.
9. Run `09_drug_sensitivity_qPCR.R` after obtaining GDSC2/oncoPredict resources and qPCR input tables.
10. Run `11_model_benchmarking_TMB_GSE50081.R`, `12_make_supportive_figures.R`, and/or `13_GSE50081_supportive_evaluation.R` for the supportive model-benchmarking, mutation-context, and GSE50081 analyses.

## Before running the scripts

- These scripts assume that input files are present in the working directory unless otherwise noted.
- The GSE68465 script includes a local filename-based GEO workflow; download the series matrix file from GEO if direct download is unavailable.
- CIBERSORTx and TIDE are external tools; their returned output files are not generated entirely inside R.
- oncoPredict/GDSC2 files are not redistributed here; obtain them from their original sources.
- RT-qPCR input tables are not bundled in this code repository because the manuscript supplementary materials provide the raw and processed RT-qPCR data.

Read the input requirements at the top of each script before running it. Some analyses also use objects created by earlier scripts in the same R session. This is a set of analysis scripts, not a single command that downloads every input and reproduces the whole study.

## License

No software license file is currently included in this repository. Third-party datasets and tools retain their own licenses and terms.
