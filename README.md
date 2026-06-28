# Analysis Code for a LUAD Hypoxia-Related Nine-Gene Signature

This repository contains analysis scripts supporting the manuscript:

**A Hypoxia-Related Nine-Gene Signature Associated With Overall Survival and Tumor Biology in Lung Adenocarcinoma: A Retrospective Observational Study**

The scripts are provided to support transparency and reproducibility of the retrospective public-cohort analyses, single-cell projection, pathway and immune-context analyses, drug-response estimation, genomic-context analyses, same-cohort model benchmarking, and RT-qPCR hypoxia-response support described in the manuscript.

## Repository Contents

| Script | Purpose |
|---|---|
| `01_download_TCGA_data.R` | Download and preprocess TCGA-LUAD RNA-seq and clinical data using TCGAbiolinks; extract MSigDB Hallmark hypoxia genes. |
| `02_differential_expression.R` | DESeq2 tumor-vs-normal differential expression analysis and hypoxia-related DEG visualization. |
| `03_LASSO_Cox_model.R` | Univariate Cox screening, LDHA/GAPDH exclusion, LASSO-Cox model construction, risk-score calculation, and TCGA survival plots. |
| `04_prognostic_evaluation.R` | Cox regression, nomogram, calibration, C-index comparison, and baseline clinical table generation. |
| `05_GSE68465_validation.R` | External public-cohort evaluation in GSE68465 using TCGA-derived coefficients and cohort-specific median stratification. |
| `06_scRNA_analysis.R` | Single-cell RNA-seq analysis and projection of the risk signature into LUAD tumor-derived cell populations. |
| `07_pathway_enrichment.R` | limma risk-group DEG analysis, GO/KEGG enrichment, Hallmark GSEA, and GSVA pathway activity analyses. |
| `08_immune_microenvironment.R` | Immune checkpoint expression, ESTIMATE, CIBERSORTx input/output handling, and TIDE input/output handling. |
| `09_drug_sensitivity_qPCR.R` | GDSC2/oncoPredict-based drug-response estimation and RT-qPCR plotting support. |
| `10_GEO_data_extraction.R` | GSE68465 GEO data extraction and probe-to-gene preprocessing. |
| `11_model_benchmarking_TMB_GSE50081.R` | Supportive analyses for model benchmarking, TCGA TMB/driver mutation context, and GSE50081 public-cohort evaluation. |
| `12_make_supportive_figures.R` | Auxiliary figure generation for benchmarking, TMB/driver mutation context, and GSE50081 supportive analyses. |
| `13_GSE50081_supportive_evaluation.R` | Standalone supportive public-cohort evaluation in GSE50081. |

Some filenames retain earlier internal development-stage names for traceability. The final manuscript text and figure legends should be used as the authoritative mapping for the submitted version.

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

Exact package versions may affect plotting details and online data downloads. Record `sessionInfo()` after running the scripts if exact computational provenance is required.

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

## Important Notes

- These scripts assume that input files are present in the working directory unless otherwise noted.
- The GSE68465 script includes a local filename-based GEO workflow; download the series matrix file from GEO if direct download is unavailable.
- CIBERSORTx and TIDE are external tools; their returned output files are not generated entirely inside R.
- oncoPredict/GDSC2 files are not redistributed here; obtain them from their original sources.
- RT-qPCR input tables are not bundled in this code repository because the manuscript supplementary materials provide the raw and processed RT-qPCR data.

## Data Availability Statement

Analysis code supporting this study is publicly available at:

https://github.com/liuwenxu44-hub/LUAD-hypoxia-nine-gene-signature

Public datasets analyzed in this study are available from TCGA/GDC and GEO under TCGA-LUAD, GSE68465, GSE50081, and GSE131907. Processed bioinformatic results and RT-qPCR raw and processed data are included in the Supplementary Material.

## License

Add a license before public deposition if all authors agree. For analysis code, the MIT License is commonly used. If no license is added, reuse rights are legally ambiguous even if the repository is public.
