# ==============================================================================
# Script Name:  code_related_to_Fig_1.R
# Purpose:      reconstructs Fig 1 plots
#               !! run reconstruct_StrokeTime_Seurat_objects_from_GSE337349_submission.R first
#
# Inputs:       Seurat object "GSE337349_StrokeTime.rds"
#               Table S1. StrokeTime sample characteristics.xlsx
#
# Dependencies: Seurat, patchwork, gridExtra, readxl, BiocParallel,
#               ggplotify, ggsci, ggthemes, cowplot, patchwork, ggsankey,
#               tidyverse, scater, batchelor, scran, dittoSeq
#
# Hardware:     > 96GB free RAM recommended
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(Seurat, patchwork, gridExtra, readxl, BiocParallel,
               ggplotify, ggsci, ggthemes, cowplot, patchwork, ggsankey,
               tidyverse, scater, batchelor, scran, dittoSeq)
ncores = ceiling(parallel::detectCores()/4)
set.seed(1234)

# set variables
file_path <- "path/to/combined_StrokeTime_seurat_object"
meta_file <- "/path/to/Table S1. StrokeTime sample characteristics.xlsx"
out_dir <- "path/to/output_directory/"

object_name <- 'GSE337349_StrokeTime.rds'

setwd(file_path)
getwd()

##### initialize functions #####
plot.sankey <- function(dfsl) {
  p <- ggplot(dfsl, aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = factor(node),
    label = paste0(node, "\n(", pct, "%)"))) +
    geom_sankey(flow.alpha = .6) +
    geom_sankey_label(
      aes(
        # label = node,
        x = stage(
          x,
          after_stat = x + .075 * if_else(
            x == 5, -1, 1
          )
        ),
        hjust = 0
      ),
      size = 8 / .pt, fill = "white"
    ) +
    ggsci::scale_fill_d3(palette = "category20") +
    guides(fill = "none") +
    ggthemes::theme_map(base_size = 11) +
    theme(
      axis.text.x = element_text(face = "bold")
    ) +
    labs(
      title = ""
    )
  return(p)
}
##### END initialize functions #####

# read and transform Seurat object
obj1 <- readRDS(object_name)
obj1
sce <- as.SingleCellExperiment(obj1)

# remove Seurat object to free up memory
rm(obj1)

# Normalization
sce <- logNormCounts(sce)

# Feature selection
feat.keep <- rownames(sce)[!grepl("mt-|Ddx3y|Tsix|Xist|Rp[ls]|Hb[ab]\\-",rownames(sce))]
dec <- modelGeneVar(sce, block=sce$Repository.Sample.ID, subset.row =  feat.keep)
hvg <- getTopHVGs(dec, n=1000)

# order batches according to cell class diversity
unique(sce$Repository.Study.ID)
order.to.merge <- as.data.frame(table(sce$Repository.Study.ID,sce$cellclass)) %>%
  mutate(sorter = ifelse(Freq<100,0,1)) %>% 
  summarize(sorter = sum(sorter), .by =  Var1) %>% 
  arrange(-sorter) %>% pull("Var1") %>% as.character()
order.to.merge
sce$Repository.Study.ID <- factor(sce$Repository.Study.ID, levels = order.to.merge)

### Run batchelor correction on sample level
sce <- batchelor::correctExperiments(sce, batch=sce$Repository.Sample.ID,
                                     subset.row=hvg, correct.all=TRUE)
sce <- scater::runUMAP(sce, dimred = 'corrected', 
                       BPPARAM = BiocParallel::MulticoreParam(ncores))

# visualize
pu <- dittoDimPlot(sce, var = "cellclass", 
                   order = c("increasing"),
                   show.others = F, 
                   show.grid.lines = F, 
                   do.label = T, 
                   labels.size = 3,
                   size = 0.01, 
                   do.raster = T) + 
  ggtitle("") + NoAxes() + theme(aspect.ratio = 1, legend.position = "none") 
pu

# visualize UMAP by time group
pu.time <- dittoDimPlot(sce, var = "cellclass", 
                        split.by = "time.group",
                        order = c("increasing"),
                        show.others = T, 
                        show.grid.lines = T, 
                        do.label = F,
                        size = 0.01, 
                        do.raster = T) + 
  ggtitle("") + theme(aspect.ratio = 1, , legend.position = "none") 
pu.time

##### make model,sex,age.group sankey chart by cell counts #######
df <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.factor, as.character) %>%
  mutate(model = case_when(grepl("contra|sham|none",model) ~ "ctl",
                           T ~ model)) %>% 
  dplyr::select(model,age.group,sex) 

dfsl <- df %>%
  make_long(model,sex,age.group)

# Calculate node counts and percentages
dagg <- data.frame()
dagg <- df %>% 
  group_by(model) %>%
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=model) %>% 
  rbind(dagg)
dagg <- df %>% 
  group_by(sex) %>%
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=sex) %>% 
  rbind(dagg)
dagg <- df %>% 
  group_by(age.group) %>% 
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=age.group) %>% 
  rbind(dagg)

# Merge percentages back into the data
dfsl <- merge(dfsl, dagg, by = "node", all.x = TRUE)

# construct graph
p1 <- plot.sankey(dfsl) + scale_x_discrete(labels= c("model","sex","age group"))
p1

# make sc.platform,model,sorted.flag sankey chart
dff <- read_excel(meta_file) %>%
  dplyr::rename(Repository.Sample.ID = sample.ID) %>% 
  dplyr::select(Repository.Sample.ID,sorted.flag)

df <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.factor, as.character) %>%
  mutate(model = case_when(grepl("contra|sham|none",model) ~ "ctl",
                           T ~ model)) %>%
  left_join(dff) %>%
  dplyr::select(model,sorted.flag,sc.platform)

dfsl <- df %>%
  make_long(model,sorted.flag,sc.platform)

# Calculate node counts and percentages
dagg <- data.frame()
dagg <- df %>% 
  group_by(model) %>%
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=model) %>% 
  rbind(dagg)
dagg <- df %>% 
  group_by(sorted.flag) %>%
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=sorted.flag) %>% 
  rbind(dagg)
dagg <- df %>% 
  group_by(sc.platform) %>% 
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=sc.platform) %>% 
  rbind(dagg)

# Merge percentages back into the data
dfsl <- merge(dfsl, dagg, by = "node", all.x = TRUE)

# construct graph
p2 <- plot.sankey(dfsl) + scale_x_discrete(labels= c("model","flow sorted","sc.platform"))
p2

## make study,model,time_point sankey chart
dfs <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.factor, as.character) %>%
  mutate(model = case_when(model == "tMCAO_PBS" ~ "tMCAO",
                           T ~ model)) %>%
  dplyr::select(model,time.group,study.short) %>%
  dplyr::rename(study=study.short, time_point=time.group)

dfsl <- dfs %>%
  make_long(study,model,time_point)

# Calculate node counts and percentages
dagg <- data.frame()
dagg <- dfs %>% 
  group_by(model) %>%
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=model) %>% 
  rbind(dagg)
dagg <- dfs %>% 
  group_by(time_point) %>%
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=time_point) %>% 
  rbind(dagg)
dagg <- dfs %>% 
  group_by(study) %>% 
  tally() %>% 
  mutate(pct = round(n / sum(n)*100,1)) %>% 
  dplyr::rename(node=study) %>% 
  rbind(dagg)

# Merge percentages back into the data
dfsl <- merge(dfsl, dagg, by = "node", all.x = TRUE)

# construct graph
p3 <- ggplot(dfsl, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = factor(node),
  label = paste0(node, " (", pct, "%)"))) +
  geom_sankey(flow.alpha = .6) +
  geom_sankey_label(
    aes(x = stage(x, after_stat = x + .075 * if_else(x == 5, -1, 1)),
        hjust = 0),
    size = 8 / .pt, fill = "white"
  ) +
  scale_fill_manual(values= dittoSeq::dittoColors()) +
  guides(fill = "none") +
  theme_sankey(base_size = 11) +
  ggthemes::theme_map(base_size = 5) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10)
  ) +
  scale_x_discrete(labels= c("study","exp. model","time point"))
p3

col.vector = c("#BBD9E5","#7AE19A","#77A5D4","#DAAFB3","#DCA05D","#75906F","#D66069","#7774D5","#D6DC66","#A047E5",
               "#95EC53","#DD66CF","#D79ED8","#DAE4B9","#74DED5","#D2E3CD","#D4E69F","#E2DB48","#6BBDD2","#E6D1C7","#C2C3DB")

p4 <- as.data.frame(colData(sce)) %>% 
  group_by(cellclass,time.group) %>% 
  summarise(N=n()) %>% 
  group_by(cellclass) %>% 
  mutate(frac = N/sum(N)) %>%
  ggplot(aes(x = frac, y = cellclass, fill = time.group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=col.vector) +
  ggtitle("Cell class representation by time point") +
  scale_x_continuous(expand = c(0, 0)) +
  cowplot::theme_cowplot() +
  theme(axis.title = element_blank(),
        legend.key = element_rect(linewidth = 0.5),
        plot.title = element_text(face = "plain", size=16, hjust = 0.5),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 14))
p4

p <- ((p3 | p1 | p2)) / (pu + pu.time + p4)  + 
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = 'bold', size = 14))
p
ggsave(file.path(out_dir,"Fig_1_out.pdf"), p, width = 20, height = 10)


