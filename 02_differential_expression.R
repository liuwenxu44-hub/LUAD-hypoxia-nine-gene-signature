# ============================================================================
# 02_differential_expression.R
# DESeq2 differential expression analysis + Figure 1 (volcano, heatmap, Venn)
# Requires output from 01_download_TCGA_data.R
# ============================================================================

library(DESeq2)
library(org.Hs.eg.db)
library(msigdbr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ggvenn)
library(dplyr)
library(tibble)

# ============================================================================
# 1. DESeq2 differential expression analysis
# ============================================================================

counts_data <- read.csv("TCGA_LUAD_STAR_counts.csv",
                         row.names = 1, check.names = FALSE)
col_data    <- read.csv("TCGA_LUAD_sample_info.csv", row.names = 1)
counts_data <- counts_data[, rownames(col_data)]

# Filter low-expression genes
keep <- rowSums(counts_data >= 10) >= (0.5 * ncol(counts_data))
counts_data <- counts_data[keep, ]

col_data$tissue_type <- relevel(factor(col_data$tissue_type), ref = "Normal")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_data),
  colData   = col_data,
  design    = ~ tissue_type
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("tissue_type", "Tumor", "Normal"))
res_df <- as.data.frame(res)

# Map ENSEMBL IDs to gene symbols
res_df <- res_df %>% rownames_to_column(var = "ENSEMBL_ID")
res_df$ENSEMBL_Clean <- gsub("\\..*", "", res_df$ENSEMBL_ID)

res_df$GeneSymbol <- mapIds(org.Hs.eg.db,
  keys      = res_df$ENSEMBL_Clean,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

res_df_clean <- res_df %>%
  filter(!is.na(GeneSymbol)) %>%
  arrange(padj) %>%
  distinct(GeneSymbol, .keep_all = TRUE)

# Load hypoxia gene set and filter DEGs
hypoxia_genes <- read.csv("hallmark_hypoxia_genes.csv")

deg_hypoxia <- res_df_clean %>%
  filter(GeneSymbol %in% hypoxia_genes$GeneSymbol,
         padj < 0.05, abs(log2FoldChange) > 0.5) %>%
  arrange(padj)

write.csv(deg_hypoxia, "TCGA_LUAD_Hypoxia_DEGs.csv", row.names = TRUE)

# ============================================================================
# 2. Volcano plot (Figure 1B)
# ============================================================================

model_genes <- c("TPI1", "GPI", "ERO1A", "DDIT4", "ANGPTL4",
                 "AKAP12", "CAVIN3", "BCAN", "IGFBP1")

res_df_plot <- res_df_clean %>%
  filter(GeneSymbol %in% hypoxia_genes$GeneSymbol) %>%
  mutate(
    Significance = case_when(
      padj < 0.05 & log2FoldChange > 0.5  ~ "Up-regulated",
      padj < 0.05 & log2FoldChange < -0.5 ~ "Down-regulated",
      TRUE ~ "Not Significant"
    ),
    Label = ifelse(GeneSymbol %in% model_genes, GeneSymbol, "")
  )

p_volcano <- ggplot(res_df_plot, aes(x = log2FoldChange,
                                      y = -log10(padj))) +
  geom_point(data = filter(res_df_plot, Significance == "Not Significant"),
             aes(color = Significance), alpha = 0.5, size = 1.5) +
  geom_point(data = filter(res_df_plot, Significance != "Not Significant"),
             aes(color = Significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Up-regulated"    = "#E64B35",
    "Down-regulated"  = "#4DBBD5",
    "Not Significant" = "grey80"
  )) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "black", alpha = 0.5) +
  geom_text_repel(aes(label = Label), size = 3,
                  box.padding = 0.8, point.padding = 0.5,
                  segment.color = "grey50", max.overlaps = Inf,
                  fontface = "italic", color = "black") +
  labs(
    title    = "Differentially Expressed Hypoxia-related Genes",
    subtitle = "TCGA-LUAD Cohort (Tumor: N=540 vs. Normal: N=59)",
    x        = expression("log"[2]*"(Fold Change)"),
    y        = expression("-log"[10]*"(adjusted P value)"),
    caption  = "Thresholds: |log2FC| > 0.5, FDR < 0.05"
  ) +
  theme_classic() +
  theme(
    legend.position  = "top",
    legend.title     = element_blank(),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 11, color = "grey30"),
    plot.caption     = element_text(hjust = 1, face = "italic", color = "grey50"),
    aspect.ratio     = 1
  )

ggsave("Figure_1B_Volcano_Plot.pdf", plot = p_volcano,
       width = 6, height = 5.5, dpi = 300)

# ============================================================================
# 3. Heatmap of top 60 hypoxia DEGs (Figure 1C)
# ============================================================================

vsd <- vst(dds, blind = FALSE)
expr_vst <- assay(vsd)

annotation_col <- data.frame(Tissue = col_data$tissue_type)
rownames(annotation_col) <- rownames(col_data)

# Shorten TCGA barcodes to 12-char patient IDs
shorten <- function(x) make.unique(substr(x, 1, 12))
rownames(annotation_col) <- shorten(rownames(annotation_col))
colnames(expr_vst)       <- shorten(colnames(expr_vst))

annotation_col <- annotation_col[colnames(expr_vst), , drop = FALSE]

# Top 30 up + top 30 down
top30_up   <- head(arrange(deg_hypoxia, desc(log2FoldChange)), 30)
top30_down <- head(arrange(deg_hypoxia, log2FoldChange), 30)
top60 <- rbind(top30_up, top30_down)

target_ensembl <- deg_hypoxia$ENSEMBL_ID
DEGs_expr <- expr_vst[target_ensembl, ]
rownames(DEGs_expr) <- deg_hypoxia$GeneSymbol

mat_top60 <- DEGs_expr[top60$GeneSymbol, ]
annotation_col <- annotation_col[order(annotation_col$Tissue), , drop = FALSE]
mat_top60 <- mat_top60[, rownames(annotation_col)]

my_breaks <- seq(-2, 2, length.out = 51)
my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

pdf("Heatmap_Top60_Final.pdf", width = 7, height = 8)
pheatmap(mat_top60,
  scale = "row", color = my_colors, breaks = my_breaks,
  annotation_col = annotation_col,
  cluster_cols = FALSE, cluster_rows = TRUE,
  show_colnames = FALSE, show_rownames = TRUE,
  fontsize_row = 7,
  main = "Top 60 Differentially Expressed Hypoxia Genes")
dev.off()

# --- Sampled heatmap (all normals + 40 random tumors) ---
set.seed(123)
n_ids <- rownames(annotation_col)[annotation_col$Tissue == "Normal"]
t_ids <- rownames(annotation_col)[annotation_col$Tissue == "Tumor"]
sampled_ids <- c(n_ids, sample(t_ids, 40))
sub_anno <- annotation_col[sampled_ids, , drop = FALSE]
sub_anno <- sub_anno[order(match(sub_anno$Tissue, c("Normal", "Tumor"))), , drop = FALSE]
mat_sub <- mat_top60[, rownames(sub_anno)]
num_normal <- sum(sub_anno$Tissue == "Normal")

anno_row <- data.frame(
  Regulation = factor(rep(c("Up-regulated", "Down-regulated"), each = 30),
                      levels = c("Up-regulated", "Down-regulated"))
)
rownames(anno_row) <- rownames(mat_sub)

ann_colors <- list(
  Tissue     = c(Normal = "#00468B", Tumor = "#ED0000"),
  Regulation = c("Up-regulated" = "#E64B35", "Down-regulated" = "#4DBBD5")
)

p_heat <- pheatmap(mat_sub,
  scale = "row", color = my_colors, breaks = my_breaks,
  annotation_col = sub_anno, annotation_row = anno_row,
  annotation_colors = ann_colors,
  cluster_cols = FALSE, cluster_rows = FALSE,
  gaps_row = 30, gaps_col = num_normal,
  border_color = NA,
  show_colnames = FALSE, show_rownames = TRUE,
  fontsize_row = 8,
  main = "Heatmap of Hypoxia Genes (Up/Down Segregation)")

ggsave("Heatmap_Ultimate_Story.pdf", plot = p_heat$gtable,
       width = 7.5, height = 8.5, dpi = 300)

# ============================================================================
# 4. Focused 39-gene heatmap (Figure 1C final)
# ============================================================================

top15_up <- deg_hypoxia %>%
  filter(!GeneSymbol %in% model_genes) %>%
  arrange(desc(log2FoldChange)) %>%
  head(15) %>% pull(GeneSymbol)

top15_down <- deg_hypoxia %>%
  filter(!GeneSymbol %in% model_genes) %>%
  arrange(log2FoldChange) %>%
  head(15) %>% pull(GeneSymbol)

focus_genes <- c(top15_up, top15_down, model_genes)

focus_df <- deg_hypoxia %>%
  filter(GeneSymbol %in% focus_genes) %>%
  arrange(desc(log2FoldChange)) %>%
  mutate(Regulation = ifelse(log2FoldChange > 0,
                              "Up in tumor", "Down in tumor"))

mat_focus <- DEGs_expr[focus_df$GeneSymbol, rownames(sub_anno)]

anno_row2 <- data.frame(Regulation = focus_df$Regulation)
rownames(anno_row2) <- focus_df$GeneSymbol
anno_row2$Regulation <- factor(anno_row2$Regulation,
                                levels = c("Up in tumor", "Down in tumor"))

num_up <- sum(focus_df$Regulation == "Up in tumor")

p39 <- pheatmap(mat_focus,
  scale = "row", color = my_colors, breaks = my_breaks,
  annotation_col = sub_anno, annotation_row = anno_row2,
  annotation_colors = ann_colors,
  cluster_cols = FALSE, cluster_rows = FALSE,
  gaps_col = num_normal, gaps_row = num_up,
  show_colnames = FALSE, show_rownames = TRUE,
  fontsize_row = 9, border_color = NA,
  main = "Expression Profiles of Key Hypoxia-related Genes")

ggsave("Figure_1C_Heatmap_39_Key_Genes_Aligned.pdf",
       plot = p39$gtable, width = 8, height = 10, dpi = 300)

# ============================================================================
# 5. Venn diagram (Figure 1D)
# ============================================================================

list_hypoxia   <- hypoxia_genes$GeneSymbol
list_TCGA_DEGs <- res_df_clean %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.5) %>% pull(GeneSymbol)

my_gene_sets <- list(
  "MSigDB Hypoxia" = list_hypoxia,
  "TCGA-LUAD DEGs" = list_TCGA_DEGs
)

p_venn <- ggvenn(my_gene_sets,
  fill_color     = c("#FDB863", "#4DBBD5"),
  fill_alpha     = 0.5,
  stroke_color   = "white",
  stroke_size    = 1,
  set_name_size  = 5,
  text_size      = 6,
  show_percentage = FALSE
) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(2, 2, 2, 2, "cm"))

ggsave("Figure_1D_Venn_Diagram_Final.pdf", p_venn,
       width = 6.5, height = 5, dpi = 300)

message("02_differential_expression.R done.")
