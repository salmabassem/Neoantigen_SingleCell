# Build_naldini_seurat.R Validation dataset from GSE185991 - from scratch. 

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



data_dir <- "/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/Naldini/"                 
meta_rds <-  "GSE185991_Full_Patient_metadata.rds"
patient_col <- "orig.ident" 


gz <- list.files(data_dir, pattern = "\\.(mtx|tsv)\\.gz$", full.names = TRUE)
prefix <- unique(sub("_(barcodes|features|matrix)\\.(tsv|mtx)\\.gz$", "", basename(gz)))
prefix <- prefix[order(prefix)]


read_one <- function(pref) {
  mtx  <- file.path(data_dir, paste0(pref, "_matrix.mtx.gz"))
  feat <- file.path(data_dir, paste0(pref, "_features.tsv.gz"))
  bc   <- file.path(data_dir, paste0(pref, "_barcodes.tsv.gz"))

  # decide which column holds the gene names (2 if available, else 1)
  k <- min(2, ncol(read.delim(feat, header = FALSE, nrows = 1)))

  mat <- ReadMtx(mtx = mtx, features = feat, cells = bc, feature.column = k)

  Mpref <- sub("^.*_", "", pref)

  # tidy names
  rownames(mat) <- make.unique(rownames(mat))
  colnames(mat) <- paste0(Mpref, "_", colnames(mat))

  seu <- CreateSeuratObject(mat, project = "Naldini", assay = "RNA", min.features = 200)
  seu$sample  <- pref
  seu$patient <- sub("^.*_(M\\d+)$", "\\1", pref)  # e.g., "M03" from "..._M03"
  seu
}

objs <- lapply(prefix, read_one)

seu <- objs[[1]]
if (length(objs) > 1) {
  for (i in 2:length(objs)) {
    seu <- merge(seu, objs[[i]])
  }
}
seu <- JoinLayers(seu)  # make sure all layers are present
colnames(seu) <- sub("-1$", "", colnames(seu))  # remove suffix "-1" if present


# -------- add patient-level metadata --------
patient_meta <- readRDS(meta_rds)   # a data.frame  
patient_meta$barcode <- rownames(patient_meta) 
patient_meta$barcode <- sub("-1$", "", patient_meta$barcode)  # remove suffix "-1" if present
seu <- AddMetaData(seu, patient_meta)

# -------- basic QC (optional) --------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

pdf("Naldini_QC.pdf", width = 12, height = 8)
VlnPlot_scCustom(seu, features = c("percent.mt", "nCount_RNA", "nFeature_RNA"))
dev.off()

# filter cells
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)
print(seu)
# -------- standard workflow --------           
### Normalize

seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = VariableFeatures(object = seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.8)  

### Extract CT45 Expression 

ct45_genes <- grep("^CT45A", rownames(seu), value = TRUE)
ct45_sum <- Matrix::colMeans(GetAssayData(seu, layer = "data")[ct45_genes, ])
seu$CT45_Sum <- ct45_sum
head(seu @meta.data)

####

pdf("Naldini_CT45Expression.pdf", width = 40, height = 40)
FeaturePlot_scCustom(seu, features = c("CT45_Sum"), num_columns = 5, pt.size = 1) + 
  scale_color_viridis_c(option = "magma") + 
  theme(text = element_text(size = 30))
dev.off()

####
# -------- save --------
saveRDS(seu, file.path(data_dir, "Naldini_merged_seurat.rds"))

###############################################################################


## CT45 Expression in Naldini 


seuDx <- subset(seu, subset = Timepoint == "DX")
pdf("/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/SCTransformPlots/Naldini_CT45_ByPatient.pdf", width = 30, height = 20)
FeaturePlot_scCustom(seuDx, features = c("CT45_Sum"), split.by = "PatientID", ncol = 4, pt.size = 1, num_columns = 4) 
dev.off()



## Project Markers

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
sets <- lapply(BMsignatures, function(s) intersect(s, rownames(seu)))


## Apply Marker Scores and UCell

set.seed(1)

expr <- GetAssayData(seu, assay = DefaultAssay(seu), slot = "data")  
scores <- UCell::ScoreSignatures_UCell(
   as.matrix(expr),
  features = sets,          # named list
  name = names(sets),       # keep your set names in meta.data
  ncores = 4
)
scores <- as.data.frame(scores)
colnames(scores) <- names(sets)
scores <- scores[ colnames(seu), , drop = FALSE ]   # align to cells
seu <- AddMetaData(seu, scores)



seu <- AddModuleScore(
  object   = seu,
  features = VanScores,
  name     = colnames(VanScores)  
)

seu <- AddModuleScore(
  object   = seu,
  features = marker_list,
  name     = name_vec   # same length as features -> one column per set
)




### Add BoneMarrow Annotations

projection_path = "/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/"
ref <- readRDS(paste0(projection_path, 'BoneMarrowMap_SymphonyReference.rds'))
ref$save_uwot_path <- paste0(projection_path, 'BoneMarrowMap_uwot_model.uwot')


batchvar <- "orig.ident"
query <- map_Query(
    exp_query = seu@assays$RNA$counts   , 
    metadata_query = seu@meta.data,
    ref_obj = ref,
    vars = batchvar
)

query <- query %>% calculate_MappingError(., reference = ref, MAD_threshold = 2.5, 
                                          threshold_by_donor = TRUE, donor_key = batchvar) # threshold mapping error on a per-sample basis.

query@meta.data %>% 
  ggplot(aes(x = mapping_error_score, fill = mapping_error_QC)) + 
  geom_histogram(bins = 200) + facet_wrap(.~get(batchvar))


query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 



#### Plots #############################################################################


seu$CT45_Group <- ifelse(seu$PatientID %in% c("PT08", "PT15"), "CT45_NonExpress", NA)
seu$CT45_Group <- ifelse(seu$PatientID %in% c("PT09", "PT12"), "CT45_Express", seu$CT45_Group)

columns <- setdiff(colnames(seu@meta.data), colnames(query@meta.data)) 
seu <- AddMetaData(seu, query@meta.data[95:100])

for (nm in Reductions(query)) {
  seu[[nm]] <- query[[nm]]                 # adds a DimReduc as-is
}

pdf("/fs/scratch/PAS2917/salmabassem/Analysis/NPM1_scRNA/SCTransformPlots/Naldini_CT45Expression_Markers.pdf", width = 20, height = 15)

seuDx <- subset(seu, subset = Timepoint == "DX")
FeaturePlot_scCustom(seuDx, features = c("CT45_Sum"), split.by = "PatientID", ncol = 4, pt.size = 1, num_columns = 3) 

seuDx <- subset(seuDx, subset = CT45_Group %in% c("CT45_Express", "CT45_NonExpress"))
Idents(seuDx) <- seuDx$main_bpe
DimPlot_scCustom(seuDx, split.by = "CT45_Group", ncol = 4, pt.size = 1, num_columns = 2) 
Idents(seuDx) <- seuDx$main_hpca
DimPlot_scCustom(seuDx,  split.by = "CT45_Group", ncol = 4, pt.size = 1, num_columns = 2) 
Idents(seuDx) <- seuDx$refined_bpe
DimPlot_scCustom(seuDx,split.by = "CT45_Group", ncol = 4, pt.size = 1, num_columns = 2) 

Idents(seuDx) <- seuDx$predicted_CellType_Broad
DimPlot_scCustom(seuDx,split.by = "CT45_Group", ncol = 4, pt.size = 1, num_columns = 2) 
DimPlot_scCustom(seuDx,split.by = "CT45_Group", reduction = "umap_projected" , ncol = 4, pt.size = 1, num_columns = 2) 


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
      seuDx,
      split.by = "CT45_Group",
      features = f,
      raster   = TRUE, num_columns = 2
      # , na_cutoff = NA   # uncomment if your feature is in meta.data and you want to suppress NA filtering warnings
    )
    print(p)
  }
}

VlnPlot(seuDx, features = feats[1:8], group.by = "predicted_CellType_Broad", pt.size = 0, split.by = "CT45_Group")
VlnPlot(seuDx, features = feats[9:16], group.by = "predicted_CellType_Broad", pt.size = 0, split.by = "CT45_Group")
VlnPlot(seuDx, features = feats[17:25], group.by = "predicted_CellType_Broad", pt.size = 0, split.by = "CT45_Group")
VlnPlot(seuDx, features = feats[26:33], group.by = "predicted_CellType_Broad", pt.size = 0, split.by = "CT45_Group")
VlnPlot(seuDx, features = feats[34:42], group.by = "predicted_CellType_Broad", pt.size = 0, split.by = "CT45_Group")


ct_col <- "predicted_CellType_Broad"
grp_col <- "CT45_Group"
sig_cols <- feats


df_long <-  subset(seuDx, subset = !is.na(predicted_CellType_Broad))@meta.data %>%
  select(all_of(c(ct_col, grp_col, sig_cols))) %>%
  pivot_longer(all_of(sig_cols), names_to = "signature", values_to = "score")

avg_mat <- df_long %>%
  group_by(.data[[ct_col]], signature) %>%
  summarize(mean = mean(score, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = signature, values_from = mean) %>%
  as.data.frame()
rownames(avg_mat) <- avg_mat[[ct_col]]; avg_mat[[ct_col]] <- NULL
pheatmap(as.matrix(avg_mat), scale = "row", clustering_method = "ward.D2",
         main = "Avg signature per cell type (row-scaled)")

Idents(seuDx) <- seuDx$CT45_Group
for (ct in unique(seuDx$predicted_CellType_Broad)) {
  print(ct)
  de <- FindMarkers( subset(seuDx, subset = predicted_CellType_Broad == ct), ident.1 = "CT45_Express", ident.2 = "CT45_NonExpress", group.by = "CT45_Group")
  p <- EnhancedVolcano::EnhancedVolcano(de, lab = rownames(de), x = 'avg_log2FC', y = 'p_val_adj',
                                    title = paste0('CT45_Express vs CT45_NonExpress in ', ct),
                                    pCutoff = 1e-3, FCcutoff = 0.5, pointSize = 2.0, labSize = 4.0)
print(p)
} 

dev.off()

###########################################################################


## Save final Seurat Object

saveRDS(seu, file.path(projection_path, "Naldini/Naldini_merged_seurat_final.rds"))
