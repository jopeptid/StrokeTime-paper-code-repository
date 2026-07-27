# ==============================================================================
# Script Name:  reconstruct_StrokeTime_Seurat_objects_from_GSE337349_submission.R
# Purpose:      reconstructs the StrokeTime Seurat objects from GSE337349 submission rds files
#
# Inputs:       StrokeTime_raw_count_matrix.rds, StrokeTime_metadata.rds
#
# Dependencies: pacman, Seurat, patchwork, tidyverse, dittoSeq
#
# Hardware:     > 24GB free RAM recommended
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(Seurat, patchwork, tidyverse, dittoSeq)

file_path <- 'path/to/StrokeTime directory containing rds files'
setwd(file_path)
getwd()

# read the data
mtx <- readRDS('GSE337349_StrokeTime_raw_count_matrix.rds')
meta <- readRDS('GSE337349_StrokeTime_metadata.rds')

# check if the data are ok; test should return TRUE
cat('Raw data matrix and metadata are consistent:',identical(colnames(mtx),rownames(meta)))

# create Seurat object
seu <- CreateSeuratObject(counts = mtx, meta.data = meta)
# add umap coordinates as reduction
umap_matrix <- as.matrix(seu@meta.data[,c('umap_1','umap_2')])
seu[["umap"]] <- CreateDimReducObject(
  embeddings = umap_matrix,
  key        = "UMAP_",
  assay      = DefaultAssay(seu)
)
# save
seu@project.name <- 'GSE337349_StrokeTime'
saveRDS(seu, paste0(seu@project.name,'.rds'))

# make a list of cellclass objects
seu.lst <- SplitObject(seu, split.by = 'cellclass')
seu.lst <- lapply(seu.lst, NormalizeData)
names(seu.lst)

# load the cellclass of interest from the Seurat object list
cellclass_of_interest <- 'BAM'
# pull Seurat obejct
obj1 <- seu.lst[[cellclass_of_interest]]
# baptize and save
obj1@project.name <- paste0('StrokeTime_',cellclass_of_interest)
saveRDS(obj1, paste0(obj1@project.name,'.rds'))

# visualize
my.col <- dittoColors()
pts <- c(1, 0.7, 0.5, 0.7)[findInterval(ncol(obj1), c(0, 5000, 20000, 100000))]
rast <- ifelse(ncol(obj1) > 100000, T, F) 
facets <- c("time.group","model")

p1 <- dittoDimPlot(obj1, 
                   var = "celltype", 
                   show.others = F, 
                   show.grid.lines = T, 
                   do.label = T, 
                   labels.size = 3,
                   color.panel = my.col, 
                   size = pts, 
                   do.raster = rast) + 
  NoLegend() +
  ggtitle(paste0(unique(obj1$cellclass), ' (n = ', ncol(obj1),')')) +
  theme(aspect.ratio = 1,
        plot.title = element_text(face = 'bold', hjust = 0.5)) 

p2 <- dittoDimPlot(obj1, 
                   var = "celltype", 
                   split.by = facets, 
                   show.others = T, 
                   show.grid.lines = T,
                   color.panel = my.col, 
                   size = pts/2, 
                   do.raster = rast) +  
  NoLegend() +
  theme(aspect.ratio = 1) + 
  ggtitle("")

p1 + p2
