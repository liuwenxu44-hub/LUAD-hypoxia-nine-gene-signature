# ============================================================================
# 07_pathway_enrichment.R
# GSEA, GSVA, GO/KEGG enrichment between high/low risk groups (Figure 6)
# Requires: expr_surv_clean and clin_sub from 03_LASSO_Cox_model.R,
#           TCGA_LUAD_RiskScore_Clinical.csv
# ============================================================================

library(limma)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(org.Hs.eg.db)
library(GSVA)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(gridExtra)

# ============================================================================
# 1. limma differential expression between High vs Low risk groups
# ============================================================================

common_samples <- intersect(colnames(expr_surv_clean), rownames(clin_sub))
expr_gsea <- expr_surv_clean[, common_samples]
clin_gsea <- clin_sub[common_samples, ]

group_list <- factor(clin_gsea$Risk_Group, levels = c("Low", "High"))
design <- model.matrix(~ 0 + group_list)
colnames(design) <- c("Low", "High")

fit <- lmFit(expr_gsea, design)
contrast.matrix <- makeContrasts(High_vs_Low = High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

DEG_all <- topTable(fit2, coef = "High_vs_Low", number = Inf)

# Generate ranked gene list for GSEA
gene_list <- DEG_all$logFC
names(gene_list) <- rownames(DEG_all)
gene_list <- sort(gene_list, decreasing = TRUE)

# ============================================================================
# 2. Volcano plot (Figure 6A)
# ============================================================================

highlight_genes <- c("ERO1A", "LDHA", "VEGFA", "SLC2A1", "HIF1A", "SNAI1",
                     "TWIST1", "VIMENTIN", "TPI1", "GPI", "DDIT4", "ANGPTL4",
                     "AKAP12", "CAVIN3", "BCAN", "IGFBP1")

plot_df <- DEG_all %>%
  mutate(
    GeneSymbol = rownames(DEG_all),
    Significance = case_when(
      P.Value < 0.05 & logFC > 0.5  ~ "Up in High Risk",
      P.Value < 0.05 & logFC < -0.5 ~ "Down in High Risk",
      TRUE ~ "Not Significant"
    ),
    Label = ifelse(GeneSymbol %in% highlight_genes, GeneSymbol, "")
  )

p_volcano <- ggplot(plot_df, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = filter(plot_df, Significance == "Not Significant"),
             color = "grey80", alpha = 0.5, size = 1.5) +
  geom_point(data = filter(plot_df, Significance != "Not Significant"),
             aes(color = Significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Up in High Risk" = "#E64B35",
                                 "Down in High Risk" = "#4DBBD5")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = Label), size = 3,
                  box.padding = 0.65, point.padding = 0.4,
                  segment.color = "grey50", max.overlaps = Inf,
                  fontface = "italic", color = "black") +
  labs(title = "DEGs between High and Low Risk Groups",
       x = expression("log"[2]*"(Fold Change)"),
       y = expression("-log"[10]*"(P-value)"),
       caption = "Thresholds: |log2FC| > 0.5, P-value < 0.05") +
  theme_classic() +
  theme(legend.position = "top", legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        plot.caption = element_text(hjust = 1, face = "italic", color = "grey50"),
        aspect.ratio = 1)

ggsave("Figure_6A_Volcano.pdf", plot = p_volcano,
       width = 8, height = 8, dpi = 300)

# ============================================================================
# 3. GO and KEGG enrichment of up-regulated genes (Figure 6B, 6C)
# ============================================================================

up_genes <- DEG_all %>%
  filter(logFC > 0.5 & P.Value < 0.05) %>% rownames()

gene_ids <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH",
                pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                readable = TRUE)

ekegg <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa",
                    pvalueCutoff = 0.05)
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

p_go   <- dotplot(ego, showCategory = 10, title = "GO (Biological Process)") +
  theme_bw()
p_kegg <- dotplot(ekegg, showCategory = 10, title = "KEGG Pathway") +
  theme_bw()

ggsave("Figure_6B_upKEGG.pdf", p_kegg, width = 8, height = 6, dpi = 300)
ggsave("Figure_6C_upGO.pdf", p_go, width = 8, height = 6, dpi = 300)

# Also for down-regulated genes
down_genes <- DEG_all %>%
  filter(logFC < -0.5 & P.Value < 0.05) %>% rownames()
gene_ids_dn <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

ego_dn <- enrichGO(gene = gene_ids_dn$ENTREZID, OrgDb = org.Hs.eg.db,
                   ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ekegg_dn <- enrichKEGG(gene = gene_ids_dn$ENTREZID, organism = "hsa",
                       pvalueCutoff = 0.05)
ekegg_dn <- setReadable(ekegg_dn, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

ggsave("sup_downGO.pdf", dotplot(ego_dn, showCategory = 10), width = 6, height = 8)
ggsave("sup_downKEGG.pdf", dotplot(ekegg_dn, showCategory = 10), width = 6, height = 8)

# ============================================================================
# 4. GSEA of Hallmark gene sets (Figure 6D)
# ============================================================================

library(clusterProfiler)
m_t2g <- read.gmt("h.all.v2026.1.Hs.symbols.gmt")

set.seed(123)
gsea_res <- GSEA(geneList = gene_list, TERM2GENE = m_t2g,
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1,
  pAdjustMethod = "BH", verbose = FALSE)

target_pathways <- c("HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS",
                     "HALLMARK_ANGIOGENESIS",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")

# Clean pathway names for display
clean_names <- c("HYPOXIA", "GLYCOLYSIS", "ANGIOGENESIS", "EMT")

gsea_clean <- gsea_res
for (i in seq_along(target_pathways)) {
  idx <- which(gsea_clean@result$ID == target_pathways[i])
  if (length(idx) > 0) {
    gsea_clean@result$ID[idx] <- clean_names[i]
    gsea_clean@result$Description[idx] <- clean_names[i]
  }
}
rownames(gsea_clean@result) <- gsea_clean@result$ID

names(gsea_clean@geneSets) <- gsub("HALLMARK_", "",
                                    names(gsea_clean@geneSets))
names(gsea_clean@geneSets) <- gsub("_", " ",
                                    names(gsea_clean@geneSets))
names(gsea_clean@geneSets) <- gsub("EPITHELIAL MESENCHYMAL TRANSITION",
                                    "EMT", names(gsea_clean@geneSets))

p_gsea <- gseaplot2(gsea_clean, geneSetID = clean_names,
  title = "Hallmark Pathways Enriched in the High-Risk Group",
  color = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"),
  pvalue_table = FALSE, base_size = 14, rel_heights = c(1.5, 0.5, 1))

ggsave("Figure_6D_GSEA.pdf", p_gsea, width = 10, height = 6, dpi = 300)

# GSEA results table
tab_data <- gsea_clean@result[clean_names, c("ID", "NES", "p.adjust")]
tab_data$NES <- round(tab_data$NES, 2)
tab_data$p.adjust <- signif(tab_data$p.adjust, 3)
colnames(tab_data) <- c("Pathway", "NES", "p.adjust")

pdf("Figure_6D_GSEA_Table.pdf", width = 5.5, height = 2)
grid.table(tab_data)
dev.off()

# Supplementary: Hallmark dotplot
p_dot <- dotplot(gsea_clean, showCategory = 10, split = ".sign",
  title = "Top 10 Enriched Hallmark Pathways", font.size = 11) +
  facet_grid(~.sign) + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "grey90"))

ggsave("Figure_S6_Hallmark_GSEA_Dotplot.pdf", p_dot,
       width = 9, height = 6.5, dpi = 300)

# ============================================================================
# 5. GSVA enrichment analysis (Figure 6E)
# ============================================================================

gsva_target <- c("HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS",
                 "HALLMARK_ANGIOGENESIS",
                 "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                 "HALLMARK_BILE_ACID_METABOLISM",
                 "HALLMARK_KRAS_SIGNALING_DN")

all_sets <- split(m_t2g[[2]], m_t2g[[1]])
target_sets <- all_sets[gsva_target]
names(target_sets) <- c("HYPOXIA", "GLYCOLYSIS", "ANGIOGENESIS", "EMT",
                         "Bile Acid Metabolism", "KRAS Signaling (Down)")

set.seed(123)
gsva_par <- gsvaParam(exprData = as.matrix(expr_gsea),
                       geneSets = target_sets, kcdf = "Gaussian")
gsva_matrix <- gsva(gsva_par, verbose = FALSE)

gsva_df <- as.data.frame(t(gsva_matrix))
gsva_df$Sample <- rownames(gsva_df)
gsva_df$Risk_Group <- clin_gsea[gsva_df$Sample, "Risk_Group"]

# Diverging bar plot of t-statistics
gsva_bar <- gsva_df %>%
  pivot_longer(cols = all_of(names(target_sets)),
               names_to = "Pathway", values_to = "GSVA_Score") %>%
  group_by(Pathway) %>%
  summarise(
    t_value = t.test(GSVA_Score[Risk_Group == "High"],
                     GSVA_Score[Risk_Group == "Low"])$statistic,
    p_value = t.test(GSVA_Score[Risk_Group == "High"],
                     GSVA_Score[Risk_Group == "Low"])$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    Enriched_Group = ifelse(t_value > 0, "High-risk", "Low-risk"),
    Significance = case_when(
      p_value < 0.001 ~ "***", p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*", TRUE ~ ""
    )
  ) %>%
  arrange(t_value)

gsva_bar$Pathway <- factor(gsva_bar$Pathway, levels = gsva_bar$Pathway)

p_gsva <- ggplot(gsva_bar, aes(x = Pathway, y = t_value,
                                fill = Enriched_Group)) +
  geom_bar(stat = "identity", width = 0.7, color = "black", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = c("Low-risk" = "#4DBBD5", "High-risk" = "#E64B35")) +
  geom_text(aes(label = Significance,
                y = ifelse(t_value > 0, t_value + 0.6, t_value - 0.6)),
            size = 5, fontface = "bold", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  theme_classic() +
  labs(x = "", y = "t-value (High-risk vs Low-risk)",
       fill = "Pathway Enriched in:") +
  theme(axis.text.y = element_text(size = 11, face = "bold"),
        legend.position = "top")

ggsave("Figure_6E_GSVA_DivergingBar.pdf", p_gsva,
       width = 7, height = 5, dpi = 300)

# GSVA violin plots for core 4 pathways
gsva_vln <- gsva_df %>%
  pivot_longer(cols = c("HYPOXIA", "GLYCOLYSIS", "ANGIOGENESIS", "EMT"),
               names_to = "Pathway", values_to = "GSVA_Score")

gsva_vln$Risk_Group <- factor(gsva_vln$Risk_Group, levels = c("Low", "High"))
gsva_vln$Pathway <- factor(gsva_vln$Pathway,
                            levels = c("HYPOXIA", "GLYCOLYSIS", "ANGIOGENESIS", "EMT"))

n_low  <- sum(gsva_df$Risk_Group == "Low")
n_high <- sum(gsva_df$Risk_Group == "High")

p_gsva_vln <- ggviolin(gsva_vln, x = "Risk_Group", y = "GSVA_Score",
  fill = "Risk_Group", palette = c("#4DBBD5", "#E64B35"),
  add = "boxplot", add.params = list(fill = "white", width = 0.15)) +
  facet_wrap(~ Pathway, scales = "free_y", ncol = 4) +
  stat_compare_means(comparisons = list(c("Low", "High")),
                     method = "wilcox.test", label = "p.signif") +
  scale_x_discrete(labels = c("Low" = paste0("Low-risk\n(n=", n_low, ")"),
                               "High" = paste0("High-risk\n(n=", n_high, ")"))) +
  theme_bw() +
  labs(y = "Pathway GSVA Score", x = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12, color = "black"),
        strip.background = element_rect(fill = "grey90", color = "black"),
        panel.grid = element_blank())

ggsave("Figure_6E_GSVA_Violin.pdf", p_gsva_vln,
       width = 10.5, height = 4.5, dpi = 300)

message("07_pathway_enrichment.R done.")
