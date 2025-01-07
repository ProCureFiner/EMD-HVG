library(optparse)
library(Seurat)
library(dplyr)

args=commandArgs(T)
project =args[1]
sample_name = args[2]
dir = args[3]
#dir="/media/disk1/SYSU/project/bin/KernelCluster_pipeline_github/test/"
#sample_name="151673"
in_dir=paste0(dir,"/input/") 
out_path=paste0(dir,"/output/") 
file_h5=paste0(in_dir,"/",sample_name,"_filtered_feature_bc_matrix.h5")
matrix.data <- Read10X_h5(file_h5) 
#seurat_object
sample_seurat <- CreateSeuratObject(counts = matrix.data, project = project,assay = "Spatial")
#add_image
img <- Seurat::Read10X_Image(image.dir = in_dir,image.name = "tissue_lowres_image.png")
Seurat::DefaultAssay(object = img) <- 'Spatial'

sample_seurat[['image']] <- img
out_file = paste0(out_path,"/", sample_name, ".rds")

#obtain GetTissueCoordinates
coldata <- GetTissueCoordinates(sample_seurat,
                                cols = c("row", "col", "tissue"),
                                scale = NULL)

sample_seurat$tissue <- coldata[rownames(sample_seurat@meta.data), "tissue"]

sample_seurat$tissue <- ifelse(sample_seurat$tissue == 1, 
                               "on_tissue", 
                               "not_on_tissue")
# We do a comparison of profiles between spots in tissue and not in tissue


#QC
tissue_qc <- sample_seurat@meta.data %>%
  select(-orig.ident) 

# Filter useful spots 
sample_seurat <- subset(sample_seurat, subset = tissue == "on_tissue")

# Continue with analysis
#给sample_seurat[["orig.ident"]] 操作给sample_seurat@meta.data增加"orig.ident"这一列信息，记录spot所属的样本ID。
sample_seurat[["orig.ident"]] <- sample_name

# Get mitochondrial genes percentage -
sample_seurat[["percent.mt"]] <- PercentageFeatureSet(sample_seurat, 
                                                      pattern = "^MT-")

# Process the data --------------------------------------------------------------------

# Exclude MT and ribosomal genes ------------------------------------------------------

#if(force_filter) {
  
  #Ribosomal and mitochondrial genes are taken out
  mt_genes <- row.names(sample_seurat)[grepl("^MT-", row.names(sample_seurat))] 
  rps_genes <- row.names(sample_seurat)[grepl("^RPS", row.names(sample_seurat))] 
  mrp_genes <- row.names(sample_seurat)[grepl("^MRP", row.names(sample_seurat))] 
  rpl_genes <- row.names(sample_seurat)[grepl("^RPL", row.names(sample_seurat))] 
  rb_genes <- c(rps_genes, mrp_genes, rpl_genes)
  

  #rm ribosomes or mitochondrion
  sample_seurat <- sample_seurat[!rownames(sample_seurat) %in% c(rb_genes, mt_genes), ]

  sample_seurat_mito <- sample_seurat[rownames(sample_seurat) %in% c(mt_genes), ]

  sample_seurat_rb <- sample_seurat[rownames(sample_seurat) %in% c(rb_genes), ]
  #Genes expressed in less that 10 spots are filtered。
  sample_seurat <- sample_seurat[rowSums(GetAssayData(sample_seurat, assay = "Spatial") > 0) > 10, ]
  
  sample_seurat$nFeature_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial") > 0)
  sample_seurat$nCount_Spatial_filt <- colSums(GetAssayData(sample_seurat, assay = "Spatial"))
  
# SCT transform normalization ---------------------------------------------------------
sample_seurat <- SCTransform(sample_seurat, 
                             assay = "Spatial", 
                             verbose = FALSE)
#DefaultAssay(sample_seurat) <- "SCT" 
# cpm normalization ---------------------------------------------------------
sample_seurat <- NormalizeData(sample_seurat, 
                                normalization.method = 'LogNormalize', 
                                scale.factor = 10000, 
                                verbose = FALSE)

# also run standard log normalization for comparison
sample_seurat <- NormalizeData(sample_seurat, verbose = FALSE, assay = "Spatial")

DefaultAssay(sample_seurat) <- "Spatial"

sample_seurat<-FindVariableFeatures(sample_seurat) 
sample_seurat <- ScaleData(sample_seurat, 
                           verbose = FALSE, 
                           features = rownames(sample_seurat)) %>%
  Seurat::RunPCA() %>%
  Seurat::RunUMAP(reduction = 'pca', dims = 1:30, verbose = FALSE)

saveRDS(sample_seurat,out_file)
