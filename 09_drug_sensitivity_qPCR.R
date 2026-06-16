# ============================================================================
# 09_drug_sensitivity_qPCR.R
# GDSC2 drug sensitivity prediction and qPCR validation (Figure 8)
# Requires: GDSC2_Expr.rds, GDSC2_Res.rds, expr_clean & clin_gsea from prior scripts,
#           qPCR input CSV files
# ============================================================================

library(oncoPredict)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(openxlsx)

# ============================================================================
# 1. Drug sensitivity prediction via oncoPredict
# ============================================================================

GDSC2_Expr <- readRDS("GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res  <- readRDS("GDSC2_Res.rds")

calcPhenotype(
  trainingExprData      = GDSC2_Expr,
  trainingPtype         = GDSC2_Res,
  testExprData          = as.matrix(expr_clean),
  batchCorrect          = "eb",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  minNumSamples         = 10,
  printOutput           = TRUE,
  removeLowVaringGenesFrom = "homogenizeData"
)

# ============================================================================
# 2. Process predicted IC50 matrix
# ============================================================================

ic50_matrix <- read.csv("calcPhenotype_Output/DrugPredictions.csv",
                         row.names = 1, check.names = FALSE)

# Build drug name lookup
drug_info <- data.frame(
  Drug_ID   = colnames(ic50_matrix),
  Drug_Name = sub("_[0-9]+$", "", colnames(ic50_matrix)),
  stringsAsFactors = FALSE
)

# Target drugs of interest
target_drugs <- c(
  "Cisplatin", "Paclitaxel", "Docetaxel", "Gemcitabine", "Pemetrexed",
  "Erlotinib", "Gefitinib", "Afatinib", "Osimertinib",
  "Trametinib", "Selumetinib", "Binimetinib",
  "Rapamycin", "Everolimus", "Dactolisib", "MK-2206",
  "Vorinostat", "Panobinostat", "Entinostat"
)

common_samples <- intersect(rownames(ic50_matrix), rownames(clin_gsea))
df <- as.data.frame(ic50_matrix[common_samples, , drop = FALSE],
                     check.names = FALSE)
df$Sample     <- rownames(df)
df$Risk_Group <- clin_gsea[common_samples, "Risk_Group"]

df_long <- df %>%
  pivot_longer(cols = -c(Sample, Risk_Group),
               names_to = "Drug_ID", values_to = "IC50") %>%
  left_join(drug_info, by = "Drug_ID") %>%
  mutate(Risk_Group = factor(Risk_Group, levels = c("Low", "High")))

# ============================================================================
# 3. Statistical testing per drug
# ============================================================================

safe_wilcox <- function(x, group) {
  keep <- !is.na(x) & !is.na(group)
  x <- x[keep]; group <- droplevels(group[keep])
  if (length(unique(group)) != 2 || any(table(group) < 2)) return(NA_real_)
  tryCatch(wilcox.test(x ~ group)$p.value, error = function(e) NA_real_)
}

drug_stats <- df_long %>%
  filter(Drug_Name %in% target_drugs) %>%
  group_by(Drug_ID, Drug_Name) %>%
  summarise(
    n_low  = sum(Risk_Group == "Low" & !is.na(IC50)),
    n_high = sum(Risk_Group == "High" & !is.na(IC50)),
    mean_low   = mean(IC50[Risk_Group == "Low"], na.rm = TRUE),
    mean_high  = mean(IC50[Risk_Group == "High"], na.rm = TRUE),
    median_low  = median(IC50[Risk_Group == "Low"], na.rm = TRUE),
    median_high = median(IC50[Risk_Group == "High"], na.rm = TRUE),
    delta_median = median_high - median_low,
    p_val = safe_wilcox(IC50, Risk_Group),
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_val, method = "BH"),
    Sensitivity_direction = case_when(
      median_high < median_low ~ "High-risk more sensitive",
      median_high > median_low ~ "High-risk less sensitive",
      TRUE ~ "No difference"
    ),
    Significance = case_when(
      is.na(p_val) ~ "NA", p_val < 0.0001 ~ "****",
      p_val < 0.001 ~ "***", p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*", TRUE ~ "ns"
    )
  ) %>%
  arrange(p_val)

# ============================================================================
# 4. Figure 8A: selected drug sensitivity boxplots
# ============================================================================

figure_drug_ids <- c("Docetaxel_1007", "Trametinib_1372",
                     "Paclitaxel_1080", "Gefitinib_1010")

plot_data <- df_long %>%
  filter(Drug_ID %in% figure_drug_ids) %>%
  mutate(Drug_Label = paste0(Drug_Name, " (", Drug_ID, ")"),
         Drug_Label = factor(Drug_Label,
           levels = c("Docetaxel (Docetaxel_1007)",
                      "Gefitinib (Gefitinib_1010)",
                      "Paclitaxel (Paclitaxel_1080)",
                      "Trametinib (Trametinib_1372)")))

p_drug <- ggplot(plot_data, aes(x = Risk_Group, y = IC50, fill = Risk_Group)) +
  geom_boxplot(alpha = 0.85, outlier.shape = 21, outlier.size = 1.0,
               width = 0.55, color = "grey25", linewidth = 0.35) +
  scale_fill_manual(values = c(Low = "#4DBBD5FF", High = "#E64B35FF")) +
  facet_wrap(~ Drug_Label, scales = "free_y", ncol = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.y.npc = "top", size = 4.5) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18))) +
  labs(x = NULL, y = "Estimated ln(IC50)") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top", legend.title = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(size = 12, color = "black"))

ggsave("Figure_8A_Target_Drugs_IC50.pdf", p_drug,
       width = 6, height = 8, device = "cairo_pdf")

# ============================================================================
# 5. Spearman correlation heatmap (Figure S8)
# ============================================================================

select_drugs <- c("Docetaxel_1007", "Paclitaxel_1080",
                  "Gefitinib_1010", "Trametinib_1372")
df_cor <- data.frame(
  Risk_Score = clin_gsea[common_samples, "Risk_Score"],
  ic50_matrix[common_samples, select_drugs]
)
names(df_cor)[-1] <- c("Docetaxel", "Paclitaxel", "Gefitinib", "Trametinib")

cor_results <- lapply(names(df_cor)[-1], function(drug) {
  test <- cor.test(df_cor$Risk_Score, df_cor[[drug]],
                   method = "spearman", exact = FALSE)
  data.frame(Drug = drug, R = test$estimate, P = test$p.value)
}) %>% bind_rows() %>%
  mutate(
    Significance = case_when(
      P < 0.001 ~ "***", P < 0.01 ~ "**", P < 0.05 ~ "*", TRUE ~ ""),
    Label = paste0(round(R, 3), "\n", Significance),
    Y_Axis = "Risk Score"
  )

p_cor <- ggplot(cor_results, aes(x = Drug, y = Y_Axis, fill = R)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Label), color = "black", size = 4.5, fontface = "bold") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                       midpoint = 0, name = "Spearman R") +
  theme_minimal() +
  labs(x = "", y = "", title = "Correlation: Risk Score vs Predicted IC50") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1,
                                    color = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15))

ggsave("Figure_S8_Correlation_Heatmap.pdf", p_cor,
       width = 7, height = 2.5, dpi = 300)

# ============================================================================
# 6. Supplementary Table S8: Drug sensitivity statistics
# ============================================================================

wb <- createWorkbook()
addWorksheet(wb, "Target_drug_stats")
writeData(wb, "Target_drug_stats", drug_stats)

header_style <- createStyle(fontSize = 11, fontColour = "white",
  fgFill = "#4F81BD", halign = "center", textDecoration = "bold")
addStyle(wb, "Target_drug_stats", header_style, rows = 1,
         cols = 1:ncol(drug_stats), gridExpand = TRUE)
freezePane(wb, "Target_drug_stats", firstActiveRow = 2)
setColWidths(wb, "Target_drug_stats", cols = 1:ncol(drug_stats), widths = "auto")

saveWorkbook(wb, "Supplementary_Table_S8_Drug_Sensitivity_GDSC2.xlsx",
             overwrite = TRUE)

# ============================================================================
# 7. qPCR validation plots (Figure 8C-D)
# ============================================================================

# --- Target gene qPCR ---
qpcr_target <- read.csv("qPCR_targets_TPI1_GPI_ERO1A_DDIT4_plot_corrected.csv")
qpcr_positive <- read.csv("qPCR_positive_genes_CA9_SLC2A1_LDHA_plot_corrected.csv")

# qPCR target genes bar plot
p_qpcr_target <- ggplot(qpcr_target, aes(x = Gene, y = Expression, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6,
           color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Expression - SD, ymax = Expression + SD),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Hypoxia" = "#E64B35")) +
  theme_classic() +
  labs(x = "", y = "Relative mRNA Expression",
       title = "qPCR Validation of Target Genes") +
  theme(axis.text.x = element_text(face = "bold.italic", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure_8C_qPCR_Target_Genes.pdf", p_qpcr_target,
       width = 8, height = 6, dpi = 300)

# qPCR positive control genes
p_qpcr_pos <- ggplot(qpcr_positive, aes(x = Gene, y = Expression, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6,
           color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = Expression - SD, ymax = Expression + SD),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = c("Control" = "#4DBBD5", "Hypoxia" = "#E64B35")) +
  theme_classic() +
  labs(x = "", y = "Relative mRNA Expression",
       title = "qPCR Validation of Positive Control Genes") +
  theme(axis.text.x = element_text(face = "bold.italic", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figure_8D_qPCR_Positive_Genes.pdf", p_qpcr_pos,
       width = 8, height = 6, dpi = 300)

message("09_drug_sensitivity_qPCR.R done.")
