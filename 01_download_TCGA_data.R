# ============================================================================
# 01_download_TCGA_data.R
# TCGA-LUAD data download and preprocessing via TCGAbiolinks
# Output: expression matrices, clinical data, Hallmark hypoxia gene list
# ============================================================================

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# --- Download TCGA-LUAD RNA-seq (STAR - Counts) ---
query_luad <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type   = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query_luad, method = "api", files.per.chunk = 20)
luad_se <- GDCprepare(query_luad)

# --- Extract count / TPM / log2(TPM+1) matrices ---
count_mat  <- assay(luad_se, "unstranded")
tpm_mat    <- assay(luad_se, "tpm_unstrand")
log2_tpm   <- log2(tpm_mat + 1)

sample_info <- as.data.frame(colData(luad_se))
gene_info   <- as.data.frame(rowData(luad_se))

# Remove list-columns that interfere with write.csv
for (col in c("treatments", "sites_of_involvement",
              "primary_site", "disease_type")) {
  sample_info[[col]] <- NULL
}

write.csv(count_mat,   "TCGA_LUAD_STAR_counts.csv")
write.csv(tpm_mat,     "TCGA_LUAD_TPM.csv")
write.csv(log2_tpm,    "TCGA_LUAD_log2TPM.csv")
write.csv(sample_info, "TCGA_LUAD_sample_info.csv", row.names = FALSE)
write.csv(gene_info,   "TCGA_LUAD_gene_annotation.csv", row.names = FALSE)

# --- Download and clean clinical data ---
clinical_luad <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

is_simple <- function(x) is.atomic(x) || is.factor(x)
clinical_luad <- clinical_luad[, sapply(clinical_luad, is_simple), drop = FALSE]

# Build overall survival variables
clinical_luad$OS.time <- ifelse(
  is.na(clinical_luad$days_to_death),
  clinical_luad$days_to_last_follow_up,
  clinical_luad$days_to_death
)
clinical_luad$OS <- ifelse(clinical_luad$vital_status == "Dead", 1, 0)

write.csv(clinical_luad, "TCGA_LUAD_clinical.csv", row.names = FALSE)

# --- Separate tumor / normal, deduplicate by patient ---
barcode       <- colnames(log2_tpm)
tumor_samples <- barcode[substr(barcode, 14, 15) == "01"]
normal_samples <- barcode[substr(barcode, 14, 15) == "11"]

tumor_expr  <- log2_tpm[, tumor_samples]
normal_expr <- log2_tpm[, normal_samples]

pid <- substr(colnames(tumor_expr), 1, 12)
tumor_expr <- tumor_expr[, !duplicated(pid)]
colnames(tumor_expr) <- pid[!duplicated(pid)]

# --- Align expression with clinical data ---
shared <- intersect(colnames(tumor_expr), clinical_luad$submitter_id)
expr_surv <- tumor_expr[, shared]
clin_surv <- clinical_luad[match(shared, clinical_luad$submitter_id), ]

write.csv(expr_surv, "TCGA_LUAD_tumor_log2TPM_for_survival.csv")
write.csv(clin_surv, "TCGA_LUAD_clinical_for_survival.csv", row.names = FALSE)

# --- MSigDB Hallmark Hypoxia gene set ---
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens")
hypoxia_genes <- unique(m_df$gene_symbol[
  m_df$gs_name == "HALLMARK_HYPOXIA"
])

write.csv(data.frame(GeneSymbol = hypoxia_genes),
          "hallmark_hypoxia_genes.csv", row.names = FALSE)

message("01_download_TCGA_data.R done.")
