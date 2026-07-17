# ==============================================================================
# Script Name:  code_related_to_Fig_S1.R
# Purpose:      processes data and generates figures 
#               !! run reconstruct_StrokeTime_Seurat_objects_from_GSE337349_submission.R first
#
# Inputs:       Seurat object "GSE337349_StrokeTime.rds"
#
# Dependencies: pacman, Seurat, edgeR, SingleCellExperiment, scuttle, randomcoloR,
#               ComplexHeatmap, dittoSeq, ggplotify, patchwork, tidyverse
#
# Hardware:     > 48GB free RAM recommended
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(Seurat, edgeR, SingleCellExperiment, scuttle, randomcoloR,
               ComplexHeatmap, dittoSeq, ggplotify, patchwork, tidyverse)
set.seed(1234)

##### initialize functions #####
ScaterLogNormCounts <- function(object = object) {
  require(SingleCellExperiment, quietly = T)
  require(scuttle, quietly = T)
  mtx = object[["RNA"]]$counts
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mtx))
  sce <- scuttle::logNormCounts(sce)
  object[["RNA"]]$data <- as(logcounts(sce), 'dgCMatrix')
  return(object)
}
##### END initialize functions #####

file_path <- "path/to/combined_StrokeTime_seurat_object"
out_dir <- "path/to/output_directory/"
object_name <- 'GSE337349_StrokeTime.rds'

setwd(file_path)
getwd()

obj1 <- readRDS(object_name)

# prepare object
obj1 <- ScaterLogNormCounts(obj1)
obj1 <- FindVariableFeatures(obj1, nfeatures=2000, assay = "RNA")

##### Supplementary figure 1A #####

# subset genes on highly variable genes
feat.keep <- VariableFeatures(obj1)
feat.keep <- feat.keep[!grepl("^mt\\-|Ddx3y|Tsix|Xist|Rp[ls]|Hb[ab]\\-", feat.keep)]
objs <- obj1[feat.keep,]

# only keep clusters with >= 50 cells by study
cells.keep <- objs@meta.data %>% 
  rownames_to_column() %>% 
  group_by(study.short,cellclass) %>%
  dplyr::filter(n()>=50) %>% 
  pull(rowname)
objs <- objs[,cells.keep]
objs@meta.data <- droplevels(objs@meta.data)
table(objs$study.short,objs$cellclass)

# aggregate expression
bulk.seu <- AggregateExpression(objs, group.by = c("cellclass","study.short"), assays = "RNA", return.seurat = T)

# use edgeR for pseudobulk normalization and conversion into cpm
y <- edgeR::DGEList(counts = bulk.seu[["RNA"]]$counts, samples=bulk.seu@meta.data)
dim(y)
group <- y$samples$cellclass
keep.exprs <- filterByExpr(y, group = group)
y <- y[keep.exprs,, keep.lib.sizes=F]
dim(y)

y <- calcNormFactors(y)
y.cpm <- edgeR::cpm(y, log = F)

dim(y.cpm)

# perform correlation analysis on normalized counts
cmtx <- cor(as.matrix(y.cpm), method = "pearson")
dim(cmtx)

study_color <- randomcoloR::distinctColorPalette(k=21)
names(study_color) <- sort(unique(y$samples$study.short))
class_color <-  dittoColors()[1:17]
names(class_color) <- levels(obj1$cellclass)
annotations <- data.frame(class = factor(y$samples$cellclass, levels = unique(y$samples$cellclass)), study = y$samples$study.short)

ha <- HeatmapAnnotation(df = annotations, col = list(class = class_color, study = study_color))
col_fun <- circlize::colorRamp2(seq(0,1, length = 9),RColorBrewer::brewer.pal(n = 9, name = "YlOrBr"))
h <- Heatmap(
  cmtx,
  km = 1,
  width = ncol(cmtx)*unit(1, "mm"), 
  height = nrow(cmtx)*unit(1, "mm"),
  col = col_fun,
  clustering_method_rows = "ward.D2",
  clustering_method_columns =  "ward.D2",
  show_column_names = FALSE,
  show_row_names = FALSE,
  top_annotation = ha,
  use_raster = F,
  heatmap_legend_param = list(ncol = 1,title = "coeff."))

h <- as.ggplot(h) + theme(plot.margin = unit(c(0, 0, 0, 0), "inches"))
h
##### END Supplementary figure 1A #####


##### Supplementary figure 1B, C, D #####

feat <- c("nCount_RNA","nFeature_RNA","percent.mt")
p.lst <- list()
for (i in seq_along(feat)) {
  my.feat <- feat[i]
  df.plot <- obj1@meta.data |> dplyr::select(my.feat,cellclass)
  p <- ggplot(df.plot, aes_string(x='cellclass', y= my.feat)) +
    geom_violin(aes(fill = cellclass), alpha = 0.8, trim = T, adjust = 1, scale = "width", linewidth = 0.1) +
    guides(fill="none") + guides(color="none") + 
    scale_fill_manual(values=dittoColors()) +
    labs(title = my.feat) +
    scale_colour_manual(values = dittoColors()) +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "italic"),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"),
      axis.text.y = element_text(size=14),
      axis.title = element_blank(),
      axis.text.x.bottom = element_text(size=10, angle=90, hjust=1, vjust=0.5), # element_blank(),
      legend.position = "none",
      axis.ticks.x = element_blank()
    )
  p.lst[[i]] <- p
}
p <- wrap_plots(p.lst, ncol = 3)
p
##### END Supplementary figure 1B, C, D #####


##### Supplementary figure 1E #####

genes.show <- c("Slc1a2","Aldoc","Tmem212","Ccdc153","Kcnj13","Folr1","Mog","Cldn11","Nrxn3","Robo2",
                "Dcn","Col1a2","Vtn","Kcnj8","Tagln","Acta2","Cldn5","Slco1a4","Sall1","Tmem119","Mrc1","Pf4",
                "Ms4a6c","Lyz2","H2-Ab1","Cd74","Cxcr2","S100a8","Cd3g","Trac","Klrb1c","Gzma","Ms4a1","Igkc")
p1 <- Seurat::DotPlot(obj1, features = genes.show, group.by = 'cellclass') +
  coord_flip() + 
  scale_color_gradient(low = "#fcfcd9", high = "firebrick3") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1

##### END Supplementary figure 1E #####


