# ============================================================================
# 03_LASSO_Cox_model.R
# Univariate Cox + LASSO Cox regression for prognostic model building (Figure 2)
# Requires: TCGA_LUAD_Hypoxia_DEGs.csv, TCGA_LUAD_tumor_log2TPM_for_survival.csv,
#           TCGA_LUAD_clinical_for_survival.csv
# ============================================================================

library(survival)
library(glmnet)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(survminer)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)

# ============================================================================
# 1. Univariate Cox regression for all hypoxia DEGs
# ============================================================================

candidate_genes <- read.csv("TCGA_LUAD_Hypoxia_DEGs.csv", row.names = 1)
candidate_genes <- rownames(candidate_genes)

expr_surv <- read.csv("TCGA_LUAD_tumor_log2TPM_for_survival.csv",
                       row.names = 1, check.names = FALSE)
clin_surv <- read.csv("TCGA_LUAD_clinical_for_survival.csv", check.names = FALSE)

# Deduplicate clinical data
clin_surv <- clin_surv[!duplicated(clin_surv[, 2]), ]
rownames(clin_surv) <- clin_surv[, 2]
clin_surv <- clin_surv[, -2]

# Map ENSEMBL to gene symbols
expr_surv <- expr_surv %>% rownames_to_column(var = "ENSEMBL_ID")
expr_surv$ENSEMBL_Clean <- gsub("\\..*", "", expr_surv$ENSEMBL_ID)
expr_surv$GeneSymbol <- mapIds(org.Hs.eg.db,
  keys = expr_surv$ENSEMBL_Clean, column = "SYMBOL",
  keytype = "ENSEMBL", multiVals = "first"
)

expr_surv_clean <- expr_surv %>%
  filter(!is.na(GeneSymbol))

expr_surv_clean$mean_expr <- rowMeans(
  dplyr::select(expr_surv_clean, -ENSEMBL_ID, -ENSEMBL_Clean, -GeneSymbol),
  na.rm = TRUE
)

expr_surv_clean <- expr_surv_clean %>%
  arrange(desc(mean_expr)) %>%
  distinct(GeneSymbol, .keep_all = TRUE)

rownames(expr_surv_clean) <- expr_surv_clean$GeneSymbol
expr_surv_clean <- dplyr::select(expr_surv_clean,
  -ENSEMBL_ID, -ENSEMBL_Clean, -GeneSymbol, -mean_expr)

# Run univariate Cox for each candidate gene
common_genes   <- intersect(candidate_genes, rownames(expr_surv_clean))
common_samples <- intersect(colnames(expr_surv_clean), rownames(clin_surv))

expr_surv_sub <- expr_surv_clean[common_genes, common_samples]
clin_surv_sub <- clin_surv[common_samples, ]
expr_surv_t   <- as.data.frame(t(expr_surv_sub))

cox_results <- data.frame(
  Gene = character(), HR = numeric(),
  HR_95CI_Lower = numeric(), HR_95CI_Upper = numeric(),
  P_value = numeric(), stringsAsFactors = FALSE
)

for (gene in common_genes) {
  fit <- coxph(Surv(clin_surv_sub$OS.time / 365, clin_surv_sub$OS) ~
               expr_surv_t[, gene])
  s <- summary(fit)
  cox_results <- rbind(cox_results, data.frame(
    Gene = gene, HR = s$coefficients[, "exp(coef)"],
    HR_95CI_Lower = s$conf.int[, "lower .95"],
    HR_95CI_Upper = s$conf.int[, "upper .95"],
    P_value = s$coefficients[, "Pr(>|z|)"]
  ))
}

write.csv(cox_results, "Univariate_Cox_All_Genes.csv", row.names = FALSE)

prog_genes <- cox_results %>%
  filter(P_value < 0.05) %>%
  arrange(P_value)

write.csv(prog_genes, "Univariate_Cox_Significant_Genes.csv", row.names = FALSE)

# ============================================================================
# 2. LASSO Cox regression for feature selection
# ============================================================================

sig_genes <- prog_genes$Gene
expr_surv_sub <- expr_surv_clean[sig_genes, common_samples]

valid_samples <- rownames(clin_surv_sub)[
  !is.na(clin_surv_sub$OS.time) & !is.na(clin_surv_sub$OS) &
    clin_surv_sub$OS.time > 0
]
expr_clean <- expr_surv_sub[, valid_samples]
clin_clean <- clin_surv_sub[valid_samples, ]

x <- as.matrix(t(expr_clean))
y <- Surv(time = clin_clean$OS.time / 365, event = clin_clean$OS)

set.seed(123)
cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1)
best_lambda <- cv_fit$lambda.min

lasso_coef <- as.numeric(coef(cv_fit, s = best_lambda))
names(lasso_coef) <- rownames(coef(cv_fit, s = best_lambda))
lasso_genes <- names(lasso_coef)[lasso_coef != 0]
lasso_vals  <- lasso_coef[lasso_coef != 0]

lasso_result_df <- data.frame(
  Gene = lasso_genes, Coefficient = lasso_vals
)
write.csv(lasso_result_df, "LASSO_Cox_Selected_Genes.csv", row.names = FALSE)

# ============================================================================
# 3. Sensitivity analysis: remove LDHA & GAPDH, re-run LASSO
# ============================================================================

sens_genes <- setdiff(prog_genes$Gene, c("LDHA", "GAPDH"))
expr_sens  <- expr_surv_clean[rownames(expr_surv_clean) %in% sens_genes, ]

clin_clean2 <- clin_surv[
  !is.na(clin_surv$OS.time) & !is.na(clin_surv$OS) & clin_surv$OS.time > 0, ]
shared2 <- intersect(colnames(expr_sens), rownames(clin_clean2))

x2 <- as.matrix(t(expr_sens[, shared2]))
y2 <- Surv(time = clin_clean2[shared2, ]$OS.time / 365,
           event = clin_clean2[shared2, ]$OS)

set.seed(123)
cv_fit2 <- cv.glmnet(x2, y2, family = "cox", alpha = 1)
lasso_coef2 <- as.numeric(coef(cv_fit2, s = cv_fit2$lambda.min))
names(lasso_coef2) <- rownames(coef(cv_fit2, s = cv_fit2$lambda.min))
keep2 <- lasso_coef2 != 0

lasso_result_df2 <- data.frame(
  Gene = names(lasso_coef2)[keep2],
  Coefficient = lasso_coef2[keep2]
)
write.csv(lasso_result_df2,
  "LASSO_Cox_Selected_Genes_noLDHA_GAPDH.csv", row.names = FALSE)

# ============================================================================
# 4. Risk score calculation and KM survival curve
# ============================================================================

expr_sub  <- expr_clean[lasso_result_df$Gene, ]
risk_score <- as.numeric(t(expr_sub) %*% lasso_result_df$Coefficient)
clin_clean$Risk_Score <- risk_score
clin_clean$Risk_Group <- factor(
  ifelse(risk_score > median(risk_score), "High", "Low"),
  levels = c("Low", "High")
)

write.csv(clin_clean, "TCGA_LUAD_RiskScore_Clinical.csv", row.names = TRUE)

# KM plot
fit_km <- survfit(Surv(OS.time / 365, OS) ~ Risk_Group, data = clin_clean)

pdf("Overall_Survival_Hypoxia_Signature_TCGA_LUAD.pdf",
    width = 14, height = 10, onefile = FALSE)
print(ggsurvplot(fit_km, data = clin_clean,
  pval = TRUE, conf.int = FALSE, risk.table = TRUE,
  palette = c("blue", "red"),
  title = "Overall Survival of Hypoxia Signature in TCGA-LUAD",
  xlab = "Time (Years)", ylab = "Survival Probability"))
dev.off()

# ============================================================================
# 5. LASSO coefficient path plot (Figure 2B)
# ============================================================================

beta_matrix <- as.matrix(cv_fit$glmnet.fit$beta)
minus_log_lambda <- -log(cv_fit$glmnet.fit$lambda)

path_df <- as.data.frame(t(beta_matrix)) %>%
  mutate(MinusLogLambda = minus_log_lambda) %>%
  pivot_longer(-MinusLogLambda, names_to = "Gene", values_to = "Coefficient")

final_genes <- lasso_result_df$Gene
path_df <- path_df %>%
  mutate(Is_Final = ifelse(Gene %in% final_genes, "Selected", "Unselected"))

label_data <- path_df %>%
  filter(Gene %in% final_genes) %>%
  group_by(Gene) %>%
  slice_max(MinusLogLambda) %>%
  ungroup()

lambda_min_x <- -log(cv_fit$lambda.min)

p_path <- ggplot(path_df, aes(x = MinusLogLambda, y = Coefficient,
                               group = Gene, color = Is_Final,
                               linewidth = Is_Final, alpha = Is_Final)) +
  geom_line() +
  scale_color_manual(values = c(Selected = "#E64B35", Unselected = "grey80")) +
  scale_linewidth_manual(values = c(Selected = 1.0, Unselected = 0.4)) +
  scale_alpha_manual(values = c(Selected = 1, Unselected = 0.4)) +
  geom_vline(xintercept = lambda_min_x, linetype = "dashed",
             color = "black", linewidth = 0.8, alpha = 0.7) +
  annotate("text", x = lambda_min_x - 0.15,
           y = max(path_df$Coefficient) * 0.9,
           label = "lambda.min", angle = 90, vjust = 1,
           size = 4.5, fontface = "italic", color = "black") +
  geom_text_repel(data = label_data, aes(label = Gene),
    size = 4, fontface = "italic", nudge_x = 0.5,
    direction = "y", hjust = 0, segment.color = "grey50") +
  theme_classic() +
  labs(title = "LASSO coefficient profiles of candidate hypoxia-related prognostic genes",
       x = expression("-Log"(lambda)), y = "Coefficients") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text  = element_text(size = 11, color = "black"),
        plot.margin = margin(15, 50, 15, 15)) +
  coord_cartesian(clip = "off")

ggsave("Fig2B_LASSO_Coefficient_Profiles.pdf", p_path,
       width = 10, height = 7, dpi = 300)

# ============================================================================
# 6. LASSO cross-validation curve (Figure 2C)
# ============================================================================

cv_df <- data.frame(
  log_lambda = log(cv_fit$lambda),
  cvm  = cv_fit$cvm,
  cvup = cv_fit$cvup,
  cvlo = cv_fit$cvlo
)
lambda_1se_x <- log(cv_fit$lambda.1se)

p_cv <- ggplot(cv_df, aes(x = log_lambda, y = cvm)) +
  geom_errorbar(aes(ymin = cvlo, ymax = cvup),
                width = 0.05, color = "grey60", alpha = 0.8) +
  geom_point(color = "#E64B35", size = 2) +
  geom_vline(xintercept = lambda_min_x, linetype = "dashed",
             color = "black", linewidth = 0.8) +
  geom_vline(xintercept = lambda_1se_x, linetype = "dashed",
             color = "#4DBBD5", linewidth = 0.8) +
  theme_classic() +
  labs(title = "LASSO Cross-Validation Curve",
       x = expression("Log"(lambda)),
       y = "Partial Likelihood Deviance") +
  annotate("text", x = lambda_min_x, y = max(cv_df$cvup) * 0.98,
           label = "lambda.min", angle = 90, vjust = -0.5,
           size = 4.5, fontface = "italic", color = "black") +
  annotate("text", x = lambda_1se_x, y = max(cv_df$cvup) * 0.98,
           label = "lambda.1se", angle = 90, vjust = -0.5,
           size = 4.5, fontface = "italic", color = "#4DBBD5") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text  = element_text(size = 11, color = "black"))

ggsave("Fig2C_LASSO_CV_Curve.pdf", p_cv, width = 7.5, height = 5, dpi = 300)

# ============================================================================
# 7. Univariate Cox forest plot for 9 model genes (Figure 2A)
# ============================================================================

model_genes <- c("TPI1", "GPI", "ERO1A", "DDIT4", "ANGPTL4",
                 "AKAP12", "CAVIN3", "BCAN", "IGFBP1")

expr_9g <- as.data.frame(t(expr_surv_clean[model_genes, ]))
shared_9g <- intersect(rownames(expr_9g), rownames(clin_clean2))
cox_9g <- cbind(clin_clean2[shared_9g, c("OS.time", "OS")],
                expr_9g[shared_9g, ])

univ_res <- lapply(model_genes, function(g) {
  f <- as.formula(paste0("Surv(OS.time/365, OS) ~ `", g, "`"))
  s <- summary(coxph(f, data = cox_9g))
  data.frame(Gene = g,
    HR = round(s$coefficients[1, "exp(coef)"], 3),
    Lower = round(s$conf.int[1, "lower .95"], 3),
    Upper = round(s$conf.int[1, "upper .95"], 3),
    P_value = s$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE)
})

univ_gene_df <- do.call(rbind, univ_res) %>%
  mutate(
    P_show = ifelse(P_value < 0.001, "P < 0.001",
                    sprintf("P = %.3f", P_value)),
    HR_CI  = paste0(HR, " (", sprintf("%.3f", Lower),
                    "-", sprintf("%.3f", Upper), ")"),
    Type   = ifelse(HR > 1, "Risk", "Protective")
  ) %>% arrange(HR)

univ_gene_df$Gene <- factor(univ_gene_df$Gene, levels = univ_gene_df$Gene)
write.csv(univ_gene_df, "Fig2A_9Genes_Univariate_Cox.csv", row.names = FALSE)

x_max <- max(univ_gene_df$Upper, na.rm = TRUE)

p_forest <- ggplot(univ_gene_df, aes(x = HR, y = Gene)) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.2, color = "black", linewidth = 0.8) +
  geom_point(aes(color = Type), shape = 15, size = 4) +
  scale_color_manual(values = c(Risk = "#E64B35", Protective = "#4DBBD5")) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.8) +
  theme_classic() +
  labs(x = "Hazard Ratio (95% CI)",
       title = "Univariate Cox regression of nine model genes") +
  geom_text(aes(x = x_max * 1.1, label = P_show), hjust = 0, size = 4) +
  geom_text(aes(x = x_max * 1.4, label = HR_CI), hjust = 0, size = 4) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 12, face = "italic", color = "black"),
        plot.title   = element_text(hjust = 0.5, size = 14)) +
  coord_cartesian(xlim = c(min(univ_gene_df$Lower) * 0.9,
                            x_max * 1.8), clip = "off")

ggsave("Fig2A_9Genes_Univariate_Forest.pdf", p_forest,
       width = 10, height = 5.5, dpi = 300)

# ============================================================================
# 8. LASSO coefficient bar chart (Figure 2D)
# ============================================================================

plot_coef <- lasso_result_df %>%
  arrange(Coefficient) %>%
  mutate(Type = "Risk Factor")

plot_coef$Gene <- factor(plot_coef$Gene, levels = plot_coef$Gene)

p_bar <- ggplot(plot_coef, aes(x = Gene, y = Coefficient, fill = Type)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Risk Factor" = "#E64B35")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  labs(title = "Coefficients of Selected Hypoxia-Related Prognostic Genes",
       x = "", y = "LASSO Regression Coefficient") +
  theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.y = element_text(size = 12, face = "bold.italic", color = "black"),
        axis.text.x = element_text(size = 11, color = "black"),
        legend.position = "none")

ggsave("Fig2D_LASSO_Coef_BarChart.pdf", p_bar,
       width = 9, height = 5, dpi = 300)

# ============================================================================
# 9. Prognostic triplot (Figure 2E)
# ============================================================================

clin_sorted <- clin_clean %>%
  arrange(Risk_Score) %>%
  mutate(Patient_Rank = row_number(),
         Status_text  = ifelse(OS == 1, "Dead", "Alive"))

split_line <- sum(clin_sorted$Risk_Group == "Low")

p1 <- ggplot(clin_sorted, aes(Patient_Rank, Risk_Score, color = Risk_Group)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(Low = "#4DBBD5", High = "#E64B35")) +
  geom_vline(xintercept = split_line + 0.5, linetype = "dashed") +
  theme_classic() +
  labs(y = "Risk Score") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), legend.position = "top",
        legend.title = element_blank(),
        plot.margin = margin(0, 10, 0, 10))

p2 <- ggplot(clin_sorted, aes(Patient_Rank, OS.time / 365, color = Status_text)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c(Alive = "#4DBBD5", Dead = "#E64B35")) +
  geom_vline(xintercept = split_line + 0.5, linetype = "dashed") +
  theme_classic() +
  labs(y = "Survival Time (Years)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), legend.position = "top",
        legend.title = element_blank(),
        plot.margin = margin(0, 10, 0, 10))

mat_heat <- t(scale(t(expr_surv_clean[model_genes, rownames(clin_sorted)])))
heat_df <- as.data.frame(mat_heat) %>%
  mutate(Gene = rownames(.)) %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "Expression") %>%
  left_join(clin_sorted %>% tibble::rownames_to_column("Sample") %>%
            dplyr::select(Sample, Patient_Rank), by = "Sample")

heat_df$Gene <- factor(heat_df$Gene, levels = rev(model_genes))

p3 <- ggplot(heat_df, aes(Patient_Rank, Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradient2(low = "navy", mid = "white",
                       high = "firebrick3", midpoint = 0) +
  geom_vline(xintercept = split_line + 0.5, linetype = "dashed") +
  theme_classic() +
  labs(x = "Patients (Increasing Risk Score)", y = "Genes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_text(face = "bold.italic", size = 10),
        legend.position = "bottom", legend.key.width = unit(1.5, "cm"))

p_triplot <- p1 / p2 / p3 + plot_layout(heights = c(1, 1, 1.5))
ggsave("Fig2E_Prognostic_Triplot.pdf", p_triplot,
       width = 7.5, height = 9, dpi = 300)

message("03_LASSO_Cox_model.R done.")
