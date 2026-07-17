# ==============================================================================
# Script Name:  standard_preprocessing_of_cellranger_count_matrices.R
# Purpose:      create and preprocess Seurat object from cellranger count matrices
#
# Inputs:       Cellranger raw_feature_bc_matrix.h5 file
#               Table S1. StrokeTime sample characteristics.xlsx
#
# Dependencies: pacman, Seurat, qs2, harmony, singleCellTK, scuttle, data.table, dittoSeq, 
#               DoubletFinder, cowplot, patchwork, gridExtra, RColorBrewer, tidyverse, readxl
#
# Hardware:     > 16GB free RAM recommended
# ==============================================================================


if (!require("pacman")) install.packages("pacman")
pacman::p_load(Seurat, harmony, singleCellTK, scuttle, data.table, dittoSeq, DoubletFinder,
               cowplot, patchwork, gridExtra, RColorBrewer, tidyverse, readxl)
ncores = parallel::detectCores()/4
set.seed(1234)

options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

##### initialize functions #####
SubsetGenesV5 <- function(seu) {
  cat("Genes before:",nrow(seu),"\n")
  mtx <- LayerData(seu, assay = "RNA", layer = "counts")
  num.cells <- Matrix::rowSums(mtx > 0)
  genes.use <- names(num.cells[which(num.cells >= 5)])
  seu <- seu[genes.use,]
  # Recalculate UMI and feature counts from the current updates matrix
  seu[["nCount_RNA"]] <- colSums(seu[["RNA"]]$counts)
  seu[["nFeature_RNA"]] <- colSums(seu[["RNA"]]$counts > 0)
  cat("Genes after:",nrow(seu),"\n")
  return(seu)
}
##### END initialize functions #####

#### set path variables ##########
file_path = '/path/to/the/working/directory/containing h5 files'
myN = 'repository.ID of the study'  # e.g. GSE267240
meta_file = "/path/to/Table S1. StrokeTime sample characteristics.xlsx"
cellranger_out = 'path/to/the/top/cellranger/directory'

setwd(file_path)
getwd()

meta <- read_excel(meta_file) %>% 
  dplyr::rename(GEO=repository.ID, GSM=sample.ID) %>%
  distinct(GSM, .keep_all = T) %>% filter(GEO == myN)

file.lst <- list.files(path = cellranger_out, pattern = "raw_feature_bc_matrix.h5", recursive = T, full.names = T)
file.lst <- file.lst[gsub(".*(GSM[0-9]+).*","\\1",file.lst) %in% meta$GSM]
file.lst

h5.lst <- lapply(file.lst, function(x) Read10X_h5(x))
names(h5.lst) <- gsub(".*(GSM[0-9]+).*","\\1",file.lst)
names(h5.lst)
# inspect the structure
sapply(h5.lst, function(x) dim(x[]))

########## 1. Quality metrics ####################
QC.dir = "QC_plots/"
if(!dir.exists(QC.dir)) dir.create(QC.dir)

graphLst <- list()
for (i in 1:length(h5.lst)) {
  fname=names(h5.lst)[i]
  mtx <- h5.lst[[i]]
  ld <- data.frame(counts=colSums(mtx),features=colSums(mtx>0)) %>%
    dplyr::filter(counts >= 30) %>%
    pivot_longer(cols = everything(),names_to = "param", values_to =  "value")
  p <- ggplot(ld, aes(x=value, color=param, fill=param)) + 
    geom_vline(aes(xintercept=200), color="red", linetype="dashed", linewidth=0.5) +
    geom_vline(aes(xintercept=300), color="blue", linetype="dashed", linewidth=0.5) +
    geom_vline(aes(xintercept=500), color="darkgreen", linetype="dashed", linewidth=0.5) +
    geom_density(alpha=.1) +
    scale_x_log10(limits = c(10,100000)) +
    theme(legend.position="none") +
    ggtitle(fname)
  graphLst[[i]] <- p
}
p3 <- cowplot::plot_grid(plotlist =  graphLst, align = 'v')
p3
ggsave(paste0(QC.dir,myN,'_Summary_DGE_histogram.pdf'), p3, 
       device = 'pdf', height = 9, width = 12)
########## END Quality metrics ####################


########## 2. run emptyDrops and get QC data ####################
sce.lst <- list()
setwd(QC.dir)
getwd()
for (i in seq_along(h5.lst)) {
  fname <- names(h5.lst[i])
  mtx <- h5.lst[[i]]
  dim(mtx)
  colnames(mtx) <- paste(fname,colnames(mtx),sep = "_")
  ## reduce size of matrix to the top 100000 droplets
  num.mol <- sort(Matrix::colSums(mtx), decreasing = T)
  tail(num.mol)
  cells.use <- names(num.mol[1:100000])
  mtx <- mtx[ ,cells.use]
  ## reduction done
  metas <- meta[meta$GSM %in% fname,]
  dtm <- data.frame(matrix(ncol=1,nrow = ncol(mtx))) %>% dplyr::select(-1) %>%
    mutate(rowname=colnames(mtx)) %>% column_to_rownames() %>% cbind(., metas) %>%
    mutate(orig.ident = GSM)
  sce <- SingleCellExperiment(list(counts=mtx),colData=dtm)
  
  # filter cells using emptyDrops ###
  sce <- runDropletQC(sce, algorithms = "emptyDrops")
  names(colData(sce))
  scef <- subsetSCECols(sce, colData = '!is.na(sce$dropletUtils_emptyDrops_fdr)')
  scef <- subsetSCECols(scef, colData = 'scef$dropletUtils_emptyDrops_fdr < 0.01')
  
  scef <- runCellQC(scef,
                    algorithms = "QCMetrics",
                    mitoGeneLocation = "rownames",
                    mitoPrefix = "^mt-"
  )
  cat(fname,"preliminary lib size",dim(scef),"\n")
  scef <- subsetSCECols(scef, colData = "sum > 200")
  cat(fname,"lib size after trimming",dim(scef),"\n")
  scef <- runCellQC(scef,
                    algorithms = c("decontX_bg"),
                    background = sce
  )
  scef$sample <- fname
  # save html report to file
  reportCellQC(scef, output_file = paste0(fname,"_cellQC.html"))
  # slim down scef
  scef <- SingleCellExperiment(list(counts = counts(scef)),
                              colData = colData(scef))
  sce.lst[[i]] <- scef
}
names(sce.lst) <- names(h5.lst)
setwd(file_path)
getwd()

# Look at some stats
sapply(sce.lst, dim)
sapply(sce.lst,  function(x) summary(x@colData$sum))
sapply(sce.lst,  function(x) summary(x@colData$mito_percent))
sapply(sce.lst,  function(x) summary(x@colData$decontX_contamination_bg))

# subset on counts, mito, decontX
sce.lst.clean = lapply(sce.lst, function(x) subsetSCECols(x, colData = c("sum > 200",
                                                                         "sum < 100000",
                                                                         "mito_percent < 20",
                                                                         'decontX_contamination_bg < 0.7'
))) 
names(sce.lst.clean) = names(h5.lst)

# check matrix size and count distribution
sapply (sce.lst.clean, dim)
sapply(sce.lst.clean, function(x) {
  hist(log10(colSums(counts(x))), breaks = 100, xlim = c(1,5)) #, main = names(x))
})

# Make Seurat object
seu.lst = lapply(sce.lst.clean, function(x) CreateSeuratObject(counts = counts(x), meta.data = as.data.frame(colData(x))))

obj1 = merge(seu.lst[[1]], seu.lst[-1])
obj1 <- JoinLayers(obj1)
obj1
# Subset genes
obj1 <- SubsetGenesV5(obj1)
# baptize the object
obj1@project.name = myN
########## END run emptyDrops and get QC data ####################


########## 3. Run DoubletFinder on processed Seurat object #################
obj.lst <- SplitObject(obj1, split.by = "orig.ident") 

df.DF = NULL
for (i in 1:length(obj.lst)) {
  # print the sample we are on
  obj <- obj.lst[[i]]
  cat("\n",names(obj.lst[i]),"\n")
  DefaultAssay(obj) <- "RNA" 
  obj <- DietSeurat(obj, assays = "RNA", layers = "counts")
  
  # Pre-process seurat object with standard seurat workflow
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = F)
  
  # Find significant PCs
  stdv <- obj[["pca"]]@stdev
  sum.stdv <- sum(obj[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  
  # finish pre-processing
  obj <- RunUMAP(obj, dims = 1:min.pc)
  obj <- FindNeighbors(object = obj, dims = 1:min.pc)              
  obj <- FindClusters(object = obj, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(obj, PCs = 1:min.pc)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  # Homotypic doublet proportion estimate
  annotations <- obj$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  # calculate expected doublet rates
  nExp.poi <- 0.00053 + (7.6*10^-6 * ncol(obj))*ncol(obj) # the formula comes from a linear fit regression on 10X doublet prediction data
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  obj <- doubletFinder(seu = obj, 
                       PCs = 1:min.pc, 
                       pK = optimal.pk,
                       nExp = nExp.poi.adj)
  obj@meta.data <- obj@meta.data %>% dplyr::select(-any_of(contains("_doublet_finder"))) %>% 
    rename_with(~ gsub("_.*","_doublet_finder",.x), .cols = any_of(starts_with(c("pANN","DF"))))
  table(obj$DF.classifications_doublet_finder)
  DimPlot(obj, group.by = "DF.classifications_doublet_finder")
  
  dd <- data.frame(DF.classifications = obj$DF.classifications_doublet_finder, DF.score = obj$pANN_doublet_finder)
  df.DF <- rbind(df.DF, dd)
}

rm(obj.lst)
identical(rownames(obj1@meta.data),rownames(df.DF))

# add DoubletFinder annotations to obj1 metadata
obj1@meta.data = cbind(obj1@meta.data,df.DF)
table(obj1$DF.classifications)

# filter out doublets
obj1 <- subset(obj1, subset = DF.classifications == "Singlet")

qs_save(obj1, paste0(obj1@project.name,".qs"))
########### END Run DoubletFinder on processed Seurat object #################


########## 4. Seurat processing and visualization ###########################
# calculate percent mito
obj1[["percent.mt"]] <- PercentageFeatureSet(obj1, pattern = "^mt-")

# create a output directory
out.dir <- paste0(obj1@project.name,"_analysis/")
if(!dir.exists(out.dir)) dir.create(out.dir)

#### run harmony on RNA ############
DefaultAssay(obj1) <- "RNA"
obj1 <- DietSeurat(obj1, assays = "RNA")
obj1[["RNA"]] <- split(obj1[["RNA"]], f = obj1$orig.ident)
obj1 <- NormalizeData(obj1)
obj1 <- FindVariableFeatures(obj1, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
obj1 <- ScaleData(obj1, vars.to.regress = "percent.mt", verbose = FALSE)
obj1 <- RunPCA(obj1, verbose = FALSE, npcs = 40)
p <- ElbowPlot(obj1, ndims = 40)
ggsave(paste0(out.dir,obj1@project.name,"_ElbowPlot.pdf"), p)
obj1 <- IntegrateLayers(object = obj1, method = HarmonyIntegration, dims = 1:30, 
                        orig.reduction = "pca", new.reduction = 'harmony',
                        assay = "RNA", verbose = F)
# Find significant PCs
stdv <- obj1[["pca"]]@stdev
sum.stdv <- sum(obj1[["pca"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                     percent.stdv[2:length(percent.stdv)]) > 0.1),
            decreasing = T)[1] + 1
pc.add <- ifelse(ncol(obj1)>50000, 15, 10)
min.pc <- min(co1, co2) + pc.add
min.pc <- min(min.pc,30)

obj1 <- RunUMAP(obj1, reduction = "harmony", dims = 1:min.pc)
obj1 <- FindNeighbors(obj1, reduction = 'harmony', dims = 1:min.pc)
obj1 <- FindClusters(obj1, resolution = 0.5, random.seed = 1,  verbose = T)
obj1 <- JoinLayers(obj1)

# save Seurat object
qs_save(obj1, paste0(obj1@project.name,".qs"))

# make plots
my.col <- dittoColors()
pts <- ifelse(ncol(obj1) > 100000, 0.7, ifelse(ncol(obj1) > 20000, 0.3, ifelse(ncol(obj1) > 5000, 0.5,1)))
rast <- ifelse(ncol(obj1) > 100000, T, F) 
splitter <- c("condition")

p1 <- dittoDimPlot(obj1, var = "seurat_clusters", 
                   show.others = F, show.grid.lines = T, do.label = T, labels.size = 3,
                   color.panel = my.col, size = pts, do.raster = rast) + NoLegend() +
  theme(aspect.ratio = 1) + ggtitle("")

p2 <- dittoDimPlot(obj1, var = "orig.ident", 
                   show.others = F, show.grid.lines = T, do.label = F, order = "randomize",
                   color.panel = my.col, size = pts, do.raster = rast) +
  theme(aspect.ratio = 1) + ggtitle("")

p3 <- dittoDimPlot(obj1, var = "seurat_clusters", split.by = splitter, 
                   show.others = F, show.grid.lines = T,
                   color.panel = my.col, size = pts, do.raster = rast) +  NoLegend() +
  theme(aspect.ratio = 1) + ggtitle("")

n.cells <- data.frame(table(obj1$seurat_clusters)) %>% dplyr::rename(cluster=Var1,No.cells=Freq)
p <- (((p1 + p2) / p3) | gridExtra::tableGrob(n.cells, rows = NULL)) + plot_layout(widths = c(6,1)) + 
  plot_annotation(title = paste0(myN," (",ncol(obj1)," cells)"),
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5)))
ggsave(paste0(out.dir, obj1@project.name,"_UMAP.pdf"), p, height = 12,  width = 12)
p <- VlnPlot(obj1, features = c("nCount_RNA","nFeature_RNA","percent.mt"), 
             cols = my.col, ncol = 1, pt.size = 0)
ggsave(paste0(out.dir, obj1@project.name,'_vln.pdf'), p, height = 8, width = 12)

# find marker genes for each cluster
all.markers <- FindAllMarkers(object = obj1, assay = "RNA", only.pos = T, min.pct = 0.2, logfc.threshold = log2(1.5), 
                              test.use = "MAST", random.seed = 13, max.cells.per.ident=5000, verbose = T)
write.csv(all.markers, file = paste0(out.dir, obj1@project.name, '_AllMarkers.csv'), row.names = F)
as.data.frame(all.markers %>% group_by(cluster) %>% top_n(10, wt = avg_log2FC)) %>% 
  dplyr::arrange(cluster, desc(avg_log2FC)) %>% write.csv(., file = paste0(out.dir, obj1@project.name, '_Top10_AllMarkers.csv'), row.names = F)
as.data.frame(all.markers %>% group_by(cluster) %>% top_n(100, wt = avg_log2FC)) %>% 
  dplyr::arrange(cluster, desc(avg_log2FC)) %>% write.csv(., file = paste0(out.dir, obj1@project.name, '_Top100_AllMarkers.csv'), row.names = F)

########## END Seurat processing and visualization ###########################


