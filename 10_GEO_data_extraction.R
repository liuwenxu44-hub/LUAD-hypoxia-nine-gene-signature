# ============================================================================
# 10_GEO_data_extraction.R
# Download and preprocess GSE68465 microarray data from GEO
# ============================================================================

library(GEOquery)

# Read the series matrix file (download from GEO:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465)
matrix_file <- "GSE68465_series_matrix.txt"
eset <- getGEO(filename = matrix_file)

expr  <- exprs(eset)
pheno <- pData(eset)

write.csv(expr,  "GSE68465_expr_probe.csv")
write.csv(pheno, "GSE68465_pheno.csv")

cat("Expression matrix dimension:", dim(expr), "\n")
cat("Phenotype table dimension:", dim(pheno), "\n")
cat("Phenotype columns:\n")
print(colnames(pheno))

# ============================================================================
# Probe-to-gene mapping and clinical data extraction
# ============================================================================

library(hgu133a.db)
library(dplyr)

pheno <- read.csv("GSE68465_pheno.csv", row.names = 1, check.names = FALSE)

clinical <- data.frame(
  sample_id      = rownames(pheno),
  geo_accession  = pheno$geo_accession,
  disease_state  = pheno$`disease_state:ch1`,
  age            = as.numeric(gsub("[^0-9\\.]", "", pheno$`age:ch1`)),
  sex            = pheno$`Sex:ch1`,
  stage          = pheno$`disease_stage:ch1`,
  vital_status   = pheno$`vital_status:ch1`,
  OS_time        = as.numeric(gsub("[^0-9\\.]", "",
                                    pheno$`months_to_last_contact_or_death:ch1`)),
  stringsAsFactors = FALSE
)

clinical$OS_status <- ifelse(clinical$vital_status == "Dead", 1,
                       ifelse(clinical$vital_status == "Alive", 0, NA))

clinical <- clinical %>%
  filter(!is.na(OS_time), !is.na(OS_status), OS_time > 0)

write.csv(clinical, "GSE68465_clinical_OS.csv", row.names = FALSE)

# Probe-to-gene symbol mapping
expr <- read.csv("GSE68465_expr_probe.csv", row.names = 1, check.names = FALSE)

probe_anno <- AnnotationDbi::select(hgu133a.db,
  keys    = rownames(expr),
  columns = "SYMBOL",
  keytype = "PROBEID"
) %>%
  filter(!is.na(SYMBOL), !duplicated(PROBEID))

expr$PROBEID <- rownames(expr)
expr_anno <- merge(probe_anno, expr, by = "PROBEID")

expr_gene <- expr_anno %>%
  dplyr::select(-PROBEID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()

rownames(expr_gene) <- expr_gene$SYMBOL
expr_gene$SYMBOL <- NULL

write.csv(expr_gene, "GSE68465_expr_gene_symbol.csv")

message("GEO data extraction and preprocessing complete.")
