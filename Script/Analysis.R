# -------------------------------
# 1. Load Required Libraries
# -------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(SingleR)
  library(celldex)
})

set.seed(123)

# -------------------------------
# 2. Load Cell Ranger Data
# -------------------------------
data_dir <- "D:/Breast_Cancer_Single_Cell/filtered_feature_bc_matrix"

raw_counts <- Read10X(data.dir = data_dir)

# Create Seurat object
sce <- CreateSeuratObject(
  counts = raw_counts,
  project = "IDC_scRNA",
  min.cells = 3,
  min.features = 200
)

# -------------------------------
# 3. Quality Control
# -------------------------------
sce[["percent.mt"]] <- PercentageFeatureSet(
  sce,
  pattern = "^MT-"
)

sce <- subset(
  sce,
  subset =
    nFeature_RNA > 300 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)

# -------------------------------
# 4. SCTransform
# -------------------------------
sce <- SCTransform(
  sce,
  vars.to.regress = "percent.mt",
  verbose = FALSE
)

DefaultAssay(sce) <- "SCT"

# -------------------------------
# 5. PCA, Clustering, UMAP
# -------------------------------
sce <- RunPCA(sce, verbose = FALSE)

sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.4)

sce <- RunUMAP(sce, dims = 1:20)

DimPlot(sce, reduction = "umap", label = TRUE) + NoLegend()

# ============================================================
# 6. AUTOMATIC CELL-TYPE ANNOTATION (SingleR)
# ============================================================

# Convert Seurat to SingleCellExperiment
sce_sce <- as.SingleCellExperiment(sce, assay = "SCT")

# Load reference atlas
ref <- HumanPrimaryCellAtlasData()

# Run SingleR
singleR_pred <- SingleR(
  test = sce_sce,
  ref = ref,
  labels = ref$label.main
)

# Add annotations to Seurat metadata
sce$SingleR_label <- singleR_pred$labels

# Visualize annotation
DimPlot(
  sce,
  reduction = "umap",
  group.by = "SingleR_label",
  label = TRUE,
  repel = TRUE
) + NoLegend()

# ============================================================
# 7. EPITHELIAL CELL EXTRACTION
# ============================================================

# Identify epithelial cells using markers + annotation
epi_cells <- WhichCells(
  sce,
  expression = EPCAM > 1 | KRT8 > 1 | KRT18 > 1
)

epi <- subset(sce, cells = epi_cells)

DimPlot(epi, reduction = "umap", label = TRUE)

# ============================================================
# 8. TUMOR vs NORMAL EPITHELIAL CLASSIFICATION
# ============================================================

# Heuristic rule:
# Normal-like: low nCount_RNA + low EMT markers
# Tumor-like: high nCount_RNA or VIM+

epi$epi_state <- ifelse(
  epi$nCount_RNA > median(epi$nCount_RNA) | FetchData(epi, "VIM") > 1,
  "Tumor_Epithelial",
  "Normal_Epithelial"
)

epi$epi_state <- factor(
  epi$epi_state,
  levels = c("Normal_Epithelial", "Tumor_Epithelial")
)

DimPlot(
  epi,
  reduction = "umap",
  group.by = "epi_state",
  cols = c("steelblue", "firebrick"),
  pt.size = 1
)

# ============================================================
# 9. DIFFERENTIAL EXPRESSION: Tumor vs Normal Epithelium
# ============================================================

Idents(epi) <- "epi_state"

tumor_vs_normal <- FindMarkers(
  epi,
  ident.1 = "Tumor_Epithelial",
  ident.2 = "Normal_Epithelial",
  assay = "SCT",
  logfc.threshold = 0.25,
  min.pct = 0.25
)

write.csv(
  tumor_vs_normal,
  file = "Tumor_vs_Normal_Epithelial_DEGs.csv"
)

# Top DEGs visualization
top_genes <- rownames(
  tumor_vs_normal %>%
    arrange(desc(avg_log2FC)) %>%
    head(20)
)

DoHeatmap(
  epi,
  features = top_genes,
  group.by = "epi_state",
  assay = "SCT"
)

# ============================================================
# 10. Save Objects
# ============================================================

saveRDS(sce, "IDC_scRNA_Annotated_SCT.rds")
saveRDS(epi, "IDC_Epithelial_Tumor_Normal.rds")

sessionInfo()
