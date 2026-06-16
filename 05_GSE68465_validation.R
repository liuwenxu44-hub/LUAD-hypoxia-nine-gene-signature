# ============================================================================
# 05_GSE68465_validation.R
# External validation of the hypoxia prognostic signature in GSE68465 (Figure 4)
# Requires: GSE68465_expr_gene_symbol.csv, GSE68465_clinical_OS.csv,
#           LASSO_Cox_Selected_Genes.csv from 03_LASSO_Cox_model.R
# ============================================================================

library(survival)
library(survminer)
library(timeROC)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)

# ============================================================================
# 1. Load and prepare GSE68465 data
# ============================================================================

geo_expr <- read.csv("GSE68465_expr_gene_symbol.csv", row.names = 1,
                      check.names = FALSE)
geo_clin <- read.csv("GSE68465_clinical_OS.csv", row.names = 1)

# Align samples
common_samples <- intersect(colnames(geo_expr), rownames(geo_clin))
geo_expr <- geo_expr[, common_samples]
geo_clin <- geo_clin[common_samples, ]

# Convert OS_time from months to years
geo_clin$OS_years <- geo_clin$OS_time / 12

# ============================================================================
# 2. Apply the LASSO model to GSE68465
# ============================================================================

lasso_res <- read.csv("LASSO_Cox_Selected_Genes.csv")
model_genes <- lasso_res$Gene

geo_model_expr <- geo_expr[model_genes, ]
coefs <- lasso_res$Coefficient[match(model_genes, lasso_res$Gene)]

geo_clin$Risk_Score <- as.numeric(t(geo_model_expr) %*% coefs)
geo_clin$Risk_Group <- factor(
  ifelse(geo_clin$Risk_Score > median(geo_clin$Risk_Score), "High", "Low"),
  levels = c("Low", "High")
)

# ============================================================================
# 3. KM survival curve (Figure 4A)
# ============================================================================

fit_geo <- survfit(Surv(OS_years, OS_status) ~ Risk_Group, data = geo_clin)

pdf("Fig4A_GSE68465_KM.pdf", width = 14, height = 10, onefile = FALSE)
print(ggsurvplot(fit_geo, data = geo_clin,
  pval = TRUE, conf.int = FALSE, risk.table = TRUE,
  palette = c("blue", "red"),
  title = "Validation in GSE68465 Cohort",
  xlab = "Time (Years)", ylab = "Survival Probability"))
dev.off()

# Cox summary
cox_geo <- coxph(Surv(OS_years, OS_status) ~ Risk_Group, data = geo_clin)
print(summary(cox_geo))

# ============================================================================
# 4. Time-dependent ROC (Figure 4B)
# ============================================================================

ROC_geo <- timeROC(T = geo_clin$OS_years,
                    delta = geo_clin$OS_status,
                    marker = geo_clin$Risk_Score,
                    cause = 1, weighting = "marginal",
                    times = c(1, 3, 5), iid = TRUE)

roc_data_geo <- data.frame(
  FPR  = as.numeric(ROC_geo$FP),
  TPR  = as.numeric(ROC_geo$TP),
  Time = rep(c("1-Year", "3-Year", "5-Year"), each = nrow(ROC_geo$TP))
)

auc_labels_geo <- c(
  sprintf("1-Year AUC = %.3f", ROC_geo$AUC[1]),
  sprintf("3-Year AUC = %.3f", ROC_geo$AUC[2]),
  sprintf("5-Year AUC = %.3f", ROC_geo$AUC[3])
)

p_roc_geo <- ggplot(roc_data_geo, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(size = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(title = "Time-dependent ROC in GSE68465",
       x = "False Positive Rate", y = "True Positive Rate") +
  scale_color_manual(values = c("1-Year" = "#E64B35", "3-Year" = "#4DBBD5",
                                 "5-Year" = "#00A087"),
                     labels = auc_labels_geo) +
  theme(legend.position = c(0.7, 0.2), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Fig4B_GSE68465_ROC.pdf", p_roc_geo, width = 14, height = 10)

# ============================================================================
# 5. Prognostic triplot (Figure 4C)
# ============================================================================

clin_sorted_geo <- geo_clin %>%
  arrange(Risk_Score) %>%
  mutate(Patient_Rank = row_number())

split_geo <- sum(clin_sorted_geo$Risk_Group == "Low")

p1_geo <- ggplot(clin_sorted_geo, aes(Patient_Rank, Risk_Score,
                                       color = Risk_Group)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = split_geo, linetype = "dashed", color = "black") +
  scale_color_manual(values = c(Low = "#4DBBD5", High = "#E64B35")) +
  theme_classic() + labs(y = "Risk Score") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top")

p2_geo <- ggplot(clin_sorted_geo, aes(Patient_Rank, OS_years,
                                       color = as.factor(OS_status))) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = split_geo, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("0" = "#4DBBD5", "1" = "#E64B35"),
                     labels = c("Alive", "Dead")) +
  theme_classic() + labs(y = "Survival Time (Years)", color = "Status") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top")

heat_expr_geo <- geo_model_expr[, rownames(clin_sorted_geo)]
heat_expr_geo_scaled <- t(scale(t(heat_expr_geo)))
heat_expr_geo_scaled[heat_expr_geo_scaled > 2]  <- 2
heat_expr_geo_scaled[heat_expr_geo_scaled < -2] <- -2

heat_df_geo <- as.data.frame(heat_expr_geo_scaled) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(clin_sorted_geo %>% rownames_to_column("Sample") %>%
            dplyr::select(Sample, Patient_Rank), by = "Sample")

p3_geo <- ggplot(heat_df_geo, aes(Patient_Rank, Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white",
                       high = "#E64B35", midpoint = 0) +
  geom_vline(xintercept = split_geo, linetype = "dashed", color = "black") +
  theme_classic() +
  labs(x = "Patients (Increasing Risk Score)", y = "Genes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(face = "italic"))

p_tri_geo <- p1_geo / p2_geo / p3_geo + plot_layout(heights = c(1, 1, 1.5))
ggsave("Fig4C_GSE68465_Prognostic_Triplot.pdf", p_tri_geo,
       width = 7.5, height = 9, dpi = 300)

# ============================================================================
# 6. Gene expression boxplot by risk group (Figure 4D)
# ============================================================================

common_genes <- intersect(model_genes, rownames(geo_expr))
expr_model <- geo_expr[common_genes, rownames(geo_clin)]

box_df <- as.data.frame(t(expr_model)) %>%
  rownames_to_column("Sample") %>%
  left_join(geo_clin %>% rownames_to_column("Sample") %>%
            dplyr::select(Sample, Risk_Group), by = "Sample") %>%
  pivot_longer(all_of(common_genes), names_to = "Gene", values_to = "Expression")

box_df$Risk_Group <- factor(box_df$Risk_Group, levels = c("Low", "High"))

p_box <- ggplot(box_df, aes(x = Gene, y = Expression, fill = Risk_Group)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.size = 1,
               outlier.alpha = 0.4, position = position_dodge(0.8)) +
  scale_fill_manual(values = c(Low = "#4DBBD5", High = "#E64B35")) +
  theme_classic() +
  labs(x = NULL, y = expression("Gene Expression (log"[2]*"TPM)"),
       fill = "Risk Group",
       title = "Model Gene Expression by Risk Group in GSE68465") +
  stat_compare_means(aes(group = Risk_Group), label = "p.signif",
                     method = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                    face = "bold.italic", size = 12),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("Fig4D_GSE68465_RiskGroup_Boxplot.pdf", p_box,
       width = 10, height = 7, dpi = 300)

message("05_GSE68465_validation.R done.")
