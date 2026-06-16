# ============================================================================
# 06_scRNA_analysis.R
# Single-cell RNA-seq analysis of hypoxia signature in LUAD TME (Figure 5)
# Requires: GSE131907 data files (RDS), LASSO_Cox_Selected_Genes.csv
# ============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ============================================================================
# 1. Load and process scRNA-seq data
# ============================================================================

anno      <- read.table("GSE131907_Lung_Cancer_cell_annotation.txt.gz",
                         header = TRUE, sep = "\t", stringsAsFactors = FALSE)
expr_raw  <- readRDS("GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
expr_norm <- readRDS("GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")

rownames(anno) <- anno$Index
common_cells <- intersect(colnames(expr_raw), rownames(anno))
expr_raw_use <- expr_raw[, common_cells]
anno_use     <- anno[common_cells, ]

sc_obj <- CreateSeuratObject(
  counts = expr_raw_use, meta.data = anno_use,
  project = "LUAD_Hypoxia", min.cells = 3, min.features = 200
)

sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
sc_obj <- subset(sc_obj,
  subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)

# Subset to primary tumor samples
sc_primary <- subset(sc_obj,
  subset = Sample_Origin %in% c("tLung", "tL/B"))
saveRDS(sc_primary, "sc_primary_tlungtLBnlung.rds")

# Downsample to max 3000 cells per cell type
Idents(sc_primary) <- sc_primary$Cell_type
sc_obj_sub <- subset(sc_primary, downsample = 3000)

# Standard Seurat processing
sc_obj_sub <- NormalizeData(sc_obj_sub, normalization.method = "LogNormalize",
                             scale.factor = 10000)
sc_obj_sub <- FindVariableFeatures(sc_obj_sub, selection.method = "vst",
                                    nfeatures = 2000)
sc_obj_sub <- ScaleData(sc_obj_sub)
sc_obj_sub <- RunPCA(sc_obj_sub, features = VariableFeatures(sc_obj_sub))
sc_obj_sub <- RunUMAP(sc_obj_sub, dims = 1:20)

# ============================================================================
# 2. Project LASSO risk score into single cells
# ============================================================================

lasso_res <- read.csv("LASSO_Cox_Selected_Genes.csv")
valid_genes <- intersect(lasso_res$Gene, rownames(sc_obj_sub))
lasso_sub   <- lasso_res[lasso_res$Gene %in% valid_genes, ]
lasso_sub   <- lasso_sub[match(valid_genes, lasso_sub$Gene), ]

sc_expr <- GetAssayData(sc_obj_sub, assay = "RNA", layer = "data")[lasso_sub$Gene, ]
sc_risk_scores <- as.numeric(crossprod(as.matrix(sc_expr), lasso_sub$Coefficient))
sc_obj_sub$LASSO_Risk_Score <- sc_risk_scores

# ============================================================================
# 3. UMAP cell annotation (Figure 5A)
# ============================================================================

p_umap_anno <- DimPlot(sc_obj_sub, reduction = "umap",
  group.by = "Cell_type", label = TRUE, label.size = 4,
  repel = TRUE, pt.size = 0.6) +
  labs(title = "Cell Type Landscape in Tumor Microenvironment",
       x = "UMAP_1", y = "UMAP_2") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.text = element_text(size = 11))

ggsave("Fig5A_Cell_Annotation_UMAP.pdf", p_umap_anno,
       width = 10, height = 6, dpi = 300)

# ============================================================================
# 4. UMAP of projected hypoxia score (Figure 5B)
# ============================================================================

p_umap_risk <- FeaturePlot(sc_obj_sub, features = "LASSO_Risk_Score",
  order = TRUE, min.cutoff = "q5", max.cutoff = "q95", pt.size = 0.6) +
  scale_color_gradientn(colors = c("navy", "white", "firebrick3")) +
  labs(title = "Projected Hypoxia Prognostic Signature Score in the TME") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("Fig5B_Hypoxia_Signature_UMAP.pdf", p_umap_risk,
       width = 10, height = 6, dpi = 300)

# ============================================================================
# 5. Violin plot of hypoxia score by cell type (Figure 5C)
# ============================================================================

vln_data <- FetchData(sc_obj_sub, vars = c("LASSO_Risk_Score", "Cell_type"))

median_scores <- vln_data %>%
  group_by(Cell_type) %>%
  summarise(Median_Score = median(LASSO_Risk_Score)) %>%
  arrange(desc(Median_Score))

vln_data$Cell_type <- factor(vln_data$Cell_type, levels = median_scores$Cell_type)

p_vln <- ggplot(vln_data, aes(x = Cell_type, y = LASSO_Risk_Score,
                               fill = Cell_type)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.15, fill = "white", color = "black",
               outlier.shape = NA) +
  theme_classic() +
  labs(title = "Hypoxia Prognostic Signature",
       x = "Cell type", y = "Projected Hypoxia Signature Score") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("Fig5C_ViolinPlot.pdf", p_vln, width = 12, height = 6, dpi = 300)

# ============================================================================
# 6. Bar plot: proportion of high-score cells (Figure 5D)
# ============================================================================

threshold <- quantile(sc_obj_sub$LASSO_Risk_Score, 0.75, na.rm = TRUE)
sc_obj_sub$Score_Group <- ifelse(sc_obj_sub$LASSO_Risk_Score >= threshold,
                                  "High score", "Others")

prop_data <- sc_obj_sub@meta.data %>%
  group_by(Cell_type, Score_Group, .drop = FALSE) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cell_type) %>%
  mutate(Total_Cells = sum(Count), Proportion = Count / Total_Cells) %>%
  filter(Score_Group == "High score") %>%
  arrange(desc(Proportion))

prop_data$Cell_type <- factor(prop_data$Cell_type, levels = prop_data$Cell_type)

p_bar <- ggplot(prop_data, aes(Cell_type, Proportion, fill = Cell_type)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, alpha = 0.8) +
  geom_text(aes(label = percent(Proportion, accuracy = 0.1)),
            vjust = -0.8, size = 4, fontface = "bold") +
  scale_y_continuous(labels = percent_format(),
                     expand = expansion(mult = c(0, 0.15))) +
  theme_classic() +
  labs(title = "Proportion of High-Score Cells across Cell Types",
       subtitle = "High-score defined as top 25% of projected hypoxia signature score",
       x = "Cell Type in TME",
       y = "Proportion of High-Score Cells (%)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 11,
                                     color = "grey40"))

ggsave("Fig5D_Barplot_Proportion.pdf", p_bar, width = 10, height = 6, dpi = 300)

# ============================================================================
# 7. DotPlot of model genes (Figure 5E)
# ============================================================================

Idents(sc_obj_sub) <- "Cell_type"
p_dot <- DotPlot(sc_obj_sub, features = valid_genes, dot.scale = 8,
                  cluster.idents = TRUE) +
  scale_color_gradientn(colors = c("navy", "white", "firebrick3")) +
  labs(title = "Expression of Hypoxia Signature Genes in the TME",
       x = "Hypoxia signature genes", y = "Cell Types in TME") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                                    face = "bold.italic", size = 10),
        axis.text.y = element_text(face = "bold", size = 11),
        plot.title  = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("Fig5E_LASSO_Genes_DotPlot.pdf", p_dot, width = 10, height = 5, dpi = 300)

message("06_scRNA_analysis.R done.")
