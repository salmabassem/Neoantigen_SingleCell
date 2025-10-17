
library(Seurat)
library(scCustomize)
library(ggplot2)
library(readxl)
library(dplyr)
library(fgsea)   
library(UCell)
library(scCustomize)
library(cowplot)
library(Seurat)
library(BoneMarrowMap)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(tidyr)
library(pheatmap)
library(Matrix)
library(SeuratExtend)
library(future)
library(readxl)
library(Azimuth)

old_plan <- plan()  # remember current plan
on.exit(plan(old_plan), add = TRUE)
plan(sequential)                                  # avoid shipping globals
options(future.globals.maxSize = 200 * 1024^3) 


load("Six scRNA Batches Corrected.Rda")
Lambo <- subset(query, subset = project == "Lambo")
Lambo <- NormalizeData(Lambo)

Lam <- subset(query, subset = project == "Lambo")

obj <- DietSeurat(Lambo, assays = "RNA", counts = TRUE, data = TRUE, scale.data = FALSE)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
obj <- SCTransform(obj)
obj <- RunPCA(obj, npcs = 30, verbose = F)
obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = F
)
obj <- FindNeighbors(obj, dims = 1:30, reduction = "integrated.dr")
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30, reduction = "integrated.dr", return.model = TRUE, spread      = 3,       # a bit more global spacing
 repulsion.strength = 0.8)


Lambo <- obj
Lambo$CT45_Group <- ifelse(Lambo$orig.ident %in% c("GSM7494274_AML8", "GSM7494263_AML2 "), NA, Lambo$CT45_Group)

##################################### Transfer BM reductions ##############################

for (nm in Reductions(Lam)) {
  Lambo[[nm]] <- Lam[[nm]]                 # adds a DimReduc as-is
}

######### Adding Azimuth Cell Type annotations #############################################

DefaultAssay(Lambo) <- "RNA"
Lambo <- JoinLayers(Lambo)
Lambo <- RunAzimuth(Lambo, reference = "pbmcref")

#DefaultAssay(Lambo) <- "SCT"

######################### Project Markers#####################################

### We have Onur clusters, VanGalen signatures, and BoneMarrowMap signatures

path  <- "/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/Onur_markers.xlsx"
VanScores <- read.csv("/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/VanGalenSignatures.csv", header = TRUE)
BMsignatures <- fgsea::gmtPathways("/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/scAML_Differentiation_Stage_MarkerGenes.gmt")


## Onur markers
marks <- excel_sheets(path)
lst <- lapply(marks, function(s) read_excel(path, sheet = s))
names(lst) <- marks
df <- as.data.frame(matrix(NA, nrow = 2000, ncol = length(names(lst))))
colnames(df) <- names(lst)
for (i in  names(lst)) {
  new <- lst[[i]]
  new <- filter(new, new$avg_log2FC > 1.0)
  df[,i] <- new$gene[1:2000]
  
}

marker_list <- lapply(as.list(df), function(x) unique(na.omit(x)))
name_vec <- gsub("[^A-Za-z0-9]+", "_", names(marker_list))

## VanGalen signatures
VanScores <- VanScores[,-12]  # remove the last column which is empty
##BMsignatures
sets <- lapply(BMsignatures, function(s) intersect(s, rownames(Lambo)))


## Apply Marker Scores and UCell
set.seed(1)
expr <- GetAssayData(Lambo, assay = DefaultAssay(Lambo), slot = "data")  
scores <- UCell::ScoreSignatures_UCell(
   as.matrix(expr),
  features = sets,          # named list
  name = names(sets),       # keep your set names in meta.data
  ncores = 4
)
scores <- as.data.frame(scores)
colnames(scores) <- names(sets)
scores <- scores[ colnames(Lambo), , drop = FALSE ]   # align to cells
Lambo <- AddMetaData(Lambo, scores)

Lambo <- AddModuleScore(
  object   = Lambo,
  features = VanScores,
  name     = colnames(VanScores)  
)

Lambo <- AddModuleScore(
  object   = Lambo,
  features = marker_list,
  name     = name_vec   # same length as features -> one column per set
)


#########################################################################

#saveRDS(Lambo,  "/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/Lambo_merged_seurat_final.rds")
Lambo <- readRDS("/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/Lambo_merged_seurat_final.rds")

####################################################Plot ###################################################################################


cols_celltype <- c(
  # Myeloid
  "Early GMP"               = "#D55E00",  # vermillion
  "Late GMP"                = "#E69F00",  # orange
  "Monocyte"                = "#f5340eff",  # brown
  "Pro-Monocyte"            = "#2258ccff",  # golden brown
  "cDC"                     = "#eedf0dff",  # teal-green
  "pDC"                     = "#f0940aff",  # green
  "EoBasoMast Precursor"    = "#B8860B",  # dark goldenrod

  # Erythroid / Meg
  "MEP"                     = "#D62728",  # red
  "Early Erythroid"         = "#FB8072",  # salmon
  "Late Erythroid"          = "#1898cbff",  # dark red
  "Megakaryocyte Precursor" = "#d21ae2ff",  # maroon

  # Stem / Progenitor
  "HSC MPP"                 = "#24b0f1ff",  # black
  "LMPP"                    = "#77d41fff",  # dark gray
  "Cycling Progenitor"      = "#F0E442",  # yellow

  # Lymphoid
  "Naive T"                 = "#56B4E9",  # sky blue
  "CD4 Memory T"            = "#0072B2",  # blue
  "CD8 Memory T"            = "#2E86AB",  # medium blue
  "NK"                      = "#17BECF",  # cyan
  "Early Lymphoid"          = "#9467BD",  # purple
  "B"                       = "#CC79A7",  # magenta
  "Pro-B"                   = "#F0027F",  # hot magenta
  "Pre-B"                   = "#E377C2",  # pink
  "Plasma Cell"             = "#8A2BE2",  # blueviolet

  # Other
  "Stromal"                 = "#009E73"   # emerald
)

pdf("/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/SCTransformPlots/Lambo_CT45Expression_Markers.pdf", width = 20, height = 15)


sub <- subset(Lambo, subset = CT45_Group %in% c("CT45_Express", "CT45_NonExpress"))
sub$patient <- sub$orig.ident


Idents(sub) <- sub$predicted_CellType_Broad

DimPlot(sub, label = TRUE)
DimPlot_scCustom(sub, reduction = "umap",   group.by = "patient",  raster = TRUE)
DimPlot_scCustom(sub, reduction = "umap",  raster = TRUE, colors_use = cols_celltype)
DimPlot_scCustom(sub, reduction ="umap", split.by = "patient", raster = TRUE, colors_use = cols_celltype)
DimPlot_scCustom(sub, reduction ="umap", split.by = "CT45_Group", raster = TRUE, colors_use = cols_celltype)

DimPlot_scCustom(sub, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3)

Idents(sub) <- sub$Malignant
DimPlot(sub, label = TRUE)
DimPlot(sub, label = TRUE, split.by = "CT45_Group")

Idents(sub) <- sub$Classified_Celltype
DimPlot(sub, label = TRUE)
DimPlot(sub, label = TRUE, split.by = "CT45_Group")

Idents(sub) <- sub$seurat_clusters
DimPlot(sub, label = TRUE)
DimPlot_scCustom(sub, reduction ="umap", split.by = "patient", raster = TRUE, colors_use = cols_celltype)
DimPlot_scCustom(sub, reduction ="umap", split.by = "CT45_Group", raster = TRUE, colors_use = cols_celltype)

Idents(sub) <- sub$predicted.celltype.l1  ## Azimuth Output 
DimPlot(sub, label = TRUE)

Idents(sub) <- sub$predicted_CellType_Broad
DimPlot_scCustom(sub, reduction ="umap", split.by = "CT45_Group", raster = TRUE,  pt.size = 1, group.by = "predicted.celltype.l1")
DimPlot_scCustom(sub, reduction ="umap", split.by = "patient", raster = TRUE,  group.by = "predicted.celltype.l1")

## BM Projection 
Idents(sub) <- sub$predicted_CellType_Broad

DimPlot_scCustom(sub, split.by = "CT45_Group",  reduction = "umap_projected", raster = TRUE, pt.size = 1, colors_use = cols_celltype)
DimPlot_scCustom(sub, split.by = "CT45_Group",  reduction = "umap", raster = TRUE, pt.size = 1, colors_use = cols_celltype)
DimPlot_scCustom(sub, split.by = "CT45_Group", group.by = "seurat_clusters",   reduction = "umap", raster = TRUE, pt.size = 1)

{
  feats <- c(
    "High_mito9","High_mito_G2_S_10","GMP_like_S_7","Progenitor_like11",
    "CD14_Mono_like2","GMP_like_G2_8","Early_Eryth_like12","LT_HSC_like4",
    "LMPP_like5","DC_like3","T_NK15","B14",
    "HSC.Prog_Normal1","GMP_Normal2","Myeloid_Normal3",
    "HSC.Prog.like_Tumor4","GMP.like_Tumor5","Myeloid.like_Tumor6",
    "HSC.like_Tumor7","Promono.like_Tumor9","Monocyte.like_Tumor10","cDC.like_Tumor11", names(sets)
  )
  for (f in feats) {
    p <- FeaturePlot_scCustom(
      sub,
      split.by = "CT45_Group",
      features = f,
      reduction = "umap",
      raster   = TRUE, num_columns = 2
      # , na_cutoff = NA   # uncomment if your feature is in meta.data and you want to suppress NA filtering warnings
    )
    print(p)
  }
}

p <- ClusterDistrBar(sub$CT45_Group, sub$predicted_CellType_Broad)
print(p)

p <- ClusterDistrBar(sub$CT45_Group, sub$seurat_clusters)
print(p)

p <- ClusterDistrBar(sub$CT45_Group, sub$Classified_Celltype)
print(p)

tmp <- subset(sub, subset = predicted_CellType_Broad %in% c("Pro-Monocyte", "Monocyte"))
p <- ClusterDistrBar(tmp$CT45_Group, tmp$predicted_CellType_Broad)
print(p)

tmp <- subset(sub, subset = predicted_CellType_Broad %in% c("HSC MPP"  , "MEP", "Megakaryocyte Precursor" , "LMPP", "Early GMP", "Late GMP"))
p <- ClusterDistrBar(tmp$CT45_Group, tmp$predicted_CellType_Broad)
print(p)

tmp <- subset(sub, subset = predicted_CellType_Broad %in% c( "MEP", "Megakaryocyte Precursor" , "Early GMP", "Late GMP"))
p <- ClusterDistrBar(tmp$CT45_Group, tmp$predicted_CellType_Broad)
print(p)

tmp <- subset(sub, subset = Classified_Celltype %in% c("GMP" , "Progenitor"   , "HSC" , "CLP" ))
p <- ClusterDistrBar(tmp$CT45_Group, tmp$Classified_Celltype)
print(p)

genes <- c("ICAM1", "ICAM2", "ICAM3", "NCAM1", "VCAN", "SOX4")
Idents(sub) <- sub$predicted_CellType_Broad

DefaultAssay(sub) <- "RNA"
cells <- WhichCells(sub, idents = c( "HSC MPP" , "LMPP" , "Early GMP",  "Late GMP" , "Pro-Monocyte", "Monocyte" ))
p <- VlnPlot2(
  sub,
  features = genes,
  group.by = "CT45_Group",
  cells = cells,
  stat.method = "wilcox.test")
print(p)

dev.off()

###############################################################################################
### Convert to anndata for pyscenic #########


as.anndata(x = Lambo, file_path = "/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/Lambo", file_name = "Lambo_all_anndata.h5ad")

#############################################################################################
