# ============================================================================
# 08_immune_microenvironment.R
# Immune microenvironment analysis: ESTIMATE, CIBERSORT, TIDE,
# immune checkpoints (Figure 7)
# Requires: expr_gsea and clin_gsea from 07_pathway_enrichment.R
# ============================================================================

library(ggplot2)
library(ggpubr)
library(estimate)
library(dplyr)
library(tidyr)

# ============================================================================
# 1. Immune checkpoint gene expression (Figure 7C)
# ============================================================================

icb_genes <- intersect(c("CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT"),
                        rownames(expr_gsea))

icb_expr <- as.data.frame(t(expr_gsea[icb_genes, ]))
icb_expr$Sample <- rownames(icb_expr)
icb_expr$Risk_Group <- factor(clin_gsea[icb_expr$Sample, "Risk_Group"],
                               levels = c("Low", "High"))

n_low  <- sum(icb_expr$Risk_Group == "Low", na.rm = TRUE)
n_high <- sum(icb_expr$Risk_Group == "High", na.rm = TRUE)

plot_df_icb <- icb_expr %>%
  pivot_longer(cols = all_of(icb_genes),
               names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = icb_genes))

p_icb <- ggviolin(plot_df_icb, x = "Risk_Group", y = "Expression",
  fill = "Risk_Group", palette = c("#4DBBD5", "#E64B35"),
  add = "boxplot", add.params = list(fill = "white", width = 0.15)) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
  stat_compare_means(comparisons = list(c("Low", "High")),
                     method = "wilcox.test", label = "p.signif") +
  scale_x_discrete(labels = c("Low" = paste0("Low-risk\n(n=", n_low, ")"),
                               "High" = paste0("High-risk\n(n=", n_high, ")"))) +
  theme_bw() +
  labs(y = "Gene Expression (log2TPM)", x = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(size = 13, face = "bold.italic"),
        panel.grid = element_blank())

ggsave("Figure_7C_Immune_Checkpoints_Expression.pdf", p_icb,
       width = 8.5, height = 6.5, dpi = 300)

# ============================================================================
# 2. ESTIMATE: Stromal/Immune scores (Figure S7A)
# ============================================================================

input_gct  <- "tcga_estimate_input.gct"
output_gct <- "tcga_estimate_output.txt"

valid_genes <- rownames(expr_gsea)[!is.na(rownames(expr_gsea)) &
                                    rownames(expr_gsea) != ""]
expr_clean_est <- expr_gsea[valid_genes, ]

gct_df <- data.frame(NAME = rownames(expr_clean_est),
                      Description = rownames(expr_clean_est))
gct_df <- cbind(gct_df, expr_clean_est)

cat("#1.2\n", file = input_gct)
cat(nrow(expr_clean_est), "\t", ncol(expr_clean_est), "\n",
    file = input_gct, append = TRUE, sep = "")
write.table(gct_df, file = input_gct, append = TRUE,
            quote = FALSE, sep = "\t", row.names = FALSE)

estimateScore(input.ds = input_gct, output.ds = output_gct,
              platform = "illumina")

est_res <- read.table(output_gct, skip = 2, header = TRUE,
                       row.names = 1, check.names = FALSE)
est_res$Description <- NULL
est_res <- as.data.frame(t(est_res))
est_res <- data.frame(lapply(est_res, as.numeric),
                       row.names = rownames(est_res))
rownames(est_res) <- gsub("\\.", "-", rownames(est_res))

common_ids <- intersect(rownames(est_res), rownames(clin_gsea))
plot_df_est <- cbind(est_res[common_ids, ],
                      Risk_Cat = clin_gsea[common_ids, "Risk_Group"])
plot_df_est$Risk_Cat <- factor(
  ifelse(grepl("high", plot_df_est$Risk_Cat, ignore.case = TRUE),
         "High", "Low"),
  levels = c("Low", "High")
)

n_low  <- sum(plot_df_est$Risk_Cat == "Low")
n_high <- sum(plot_df_est$Risk_Cat == "High")

plot_df_long <- plot_df_est %>%
  dplyr::select(StromalScore, ImmuneScore, ESTIMATEScore, Risk_Cat) %>%
  pivot_longer(cols = 1:3, names_to = "Score_Type", values_to = "Score_Value") %>%
  mutate(Score_Type = factor(Score_Type,
         levels = c("StromalScore", "ImmuneScore", "ESTIMATEScore")))

p_estimate <- ggviolin(plot_df_long, x = "Risk_Cat", y = "Score_Value",
  fill = "Risk_Cat", palette = c("#4DBBD5", "#E64B35"),
  add = "boxplot", add.params = list(fill = "white", width = 0.15)) +
  facet_wrap(~ Score_Type, scales = "free_y", ncol = 3) +
  stat_compare_means(comparisons = list(c("Low", "High")),
                     method = "wilcox.test", label = "p.signif") +
  scale_x_discrete(labels = c("Low" = paste0("Low-risk\n(n=", n_low, ")"),
                               "High" = paste0("High-risk\n(n=", n_high, ")"))) +
  theme_bw() +
  labs(y = "ESTIMATE Score", x = "") +
  theme(legend.position = "none", panel.grid = element_blank(),
        axis.text.x = element_text(size = 11, face = "bold", color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black"))

ggsave("Figure_S7A_ESTIMATE_Score.pdf", p_estimate,
       width = 8.5, height = 6.5, dpi = 300)

# ============================================================================
# 3. CIBERSORT immune cell infiltration (Figure 7B)
# ============================================================================

# Export expression matrix for CIBERSORTx web tool
write.table(as.data.frame(expr_gsea), "CIBERSORTx_mixture.txt",
            sep = "\t", quote = FALSE, col.names = NA)

# Read CIBERSORTx results (run externally)
ciber_res <- read.csv("CIBERSORTx_Job1_Results.csv",
                       row.names = 1, check.names = FALSE)
ciber_res <- ciber_res[, 1:22]

common_ids <- intersect(rownames(ciber_res), rownames(clin_gsea))
plot_df_ciber <- cbind(ciber_res[common_ids, ],
                        Risk_Group = clin_gsea[common_ids, "Risk_Group"])
plot_df_ciber$Risk_Group <- factor(plot_df_ciber$Risk_Group,
                                    levels = c("Low", "High"))

# Filter to significant cell types only
ciber_long <- plot_df_ciber %>%
  pivot_longer(cols = 1:22, names_to = "Cell_Type", values_to = "Fraction")

sig_cells <- ciber_long %>%
  group_by(Cell_Type) %>%
  summarise(p_val = wilcox.test(Fraction ~ Risk_Group)$p.value,
            .groups = "drop") %>%
  filter(p_val < 0.05) %>%
  pull(Cell_Type)

plot_df_main <- ciber_long %>% filter(Cell_Type %in% sig_cells)

p_cibersort <- ggplot(plot_df_main,
                       aes(x = Cell_Type, y = Fraction, fill = Risk_Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.7,
               position = position_dodge(0.8), alpha = 0.8,
               color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c(Low = "#4DBBD5", High = "#E64B35")) +
  stat_compare_means(aes(group = Risk_Group), label = "p.signif",
                     method = "wilcox.test", size = 5) +
  theme_classic() +
  labs(x = "", y = "Relative Fraction of Immune Cells",
       fill = "Risk Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        legend.position = "top")

ggsave("Figure_7B_CIBERSORT_MainFigure.pdf", p_cibersort,
       width = 8, height = 6, dpi = 300)

# Full CIBERSORT landscape (Supplementary)
ggsave("Figure_S7B_CIBERSORT_Landscape.pdf", p_cibersort,
       width = 10, height = 6.5, dpi = 300)

# ============================================================================
# 4. TIDE analysis (Figure 7D, 7E)
# ============================================================================

# Center expression for TIDE input
expr_tide <- t(scale(t(expr_gsea), center = TRUE, scale = FALSE))
write.table(as.data.frame(expr_tide), "TIDE_input_centered.txt",
            sep = "\t", quote = FALSE, col.names = NA)

# Read TIDE results (run externally at http://tide.dfci.harvard.edu/)
tide_res <- read.table("TIDEx_results.csv", sep = "\t", header = TRUE,
                        row.names = 1, check.names = FALSE,
                        stringsAsFactors = FALSE)
rownames(tide_res) <- gsub("\\.", "-", rownames(tide_res))
rownames(tide_res) <- sub("^X", "", rownames(tide_res))

common_ids <- intersect(rownames(tide_res), rownames(clin_gsea))
plot_df_tide <- cbind(tide_res[common_ids, ],
                       Risk_Group = clin_gsea[common_ids, "Risk_Group"])
plot_df_tide$Risk_Group <- factor(plot_df_tide$Risk_Group,
                                   levels = c("Low", "High"))

n_low  <- sum(plot_df_tide$Risk_Group == "Low")
n_high <- sum(plot_df_tide$Risk_Group == "High")

# TIDE score violin
p_tide_score <- ggviolin(plot_df_tide, x = "Risk_Group", y = "TIDE",
  fill = "Risk_Group", palette = c("#4DBBD5", "#E64B35"),
  add = "boxplot", add.params = list(fill = "white", width = 0.15)) +
  stat_compare_means(comparisons = list(c("Low", "High")),
                     method = "wilcox.test", label = "p.signif") +
  scale_x_discrete(labels = c("Low" = paste0("Low-risk\n(n=", n_low, ")"),
                               "High" = paste0("High-risk\n(n=", n_high, ")"))) +
  theme_bw() +
  labs(y = "TIDE Score", x = "",
       title = "Tumor Immune Dysfunction and Exclusion") +
  theme(plot.title = element_text(size = 13, hjust = 0.5),
        legend.position = "none", panel.grid = element_blank())

ggsave("Figure_7D_TIDE_Score.pdf", p_tide_score,
       width = 6, height = 8, dpi = 300)

# Responder bar plot
plot_df_tide$Responder <- factor(
  ifelse(plot_df_tide$Responder %in% c("True", "TRUE", "1"),
         "Responder", "Non-Responder"),
  levels = c("Responder", "Non-Responder")
)

# Fisher's exact test
fisher_res <- fisher.test(table(plot_df_tide$Risk_Group,
                                 plot_df_tide$Responder))
p_label <- ifelse(fisher_res$p.value < 0.001,
                  "P < 0.001",
                  paste0("P = ", round(fisher_res$p.value, 3)))

responder_summary <- plot_df_tide %>%
  group_by(Risk_Group, Responder) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Risk_Group) %>%
  mutate(Total = sum(Count), Percent = Count / Total,
         Label_Text = paste0(sprintf("%.1f", Percent * 100), "%"))

p_responder <- ggplot(responder_summary,
                       aes(x = Risk_Group, y = Percent, fill = Responder)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.4) +
  geom_text(aes(label = Label_Text),
            position = position_stack(vjust = 0.5),
            size = 4, fontface = "bold", color = "white") +
  scale_fill_manual(values = c(Responder = "#374E8CFF",
                                "Non-Responder" = "#B05923FF")) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0, 0.1, 0)) +
  scale_x_discrete(labels = c("Low" = paste0("Low-risk\n(n=", n_low, ")"),
                               "High" = paste0("High-risk\n(n=", n_high, ")"))) +
  theme_classic() +
  labs(x = "", y = "Proportion of Patients", fill = "ICB Response",
       title = paste0("Predicted Clinical Response (Fisher's Exact, ", p_label, ")")) +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(size = 12, color = "black"),
        legend.position = "top")

ggsave("Figure_7E_Responder_Barplot.pdf", p_responder,
       width = 5.5, height = 7, dpi = 300)

message("08_immune_microenvironment.R done.")
