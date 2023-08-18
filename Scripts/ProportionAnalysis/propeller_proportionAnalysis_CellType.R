## spackle celltype proportions 

library(Seurat)
library(speckle)
library(dplyr)
library(ggplot2)

kelly.colours <- c("gray95", "gray13", "gold2", "plum4", "darkorange1", "lightskyblue2", "firebrick", "burlywood3", "gray51", "springgreen4", "lightpink2", "deepskyblue4", "lightsalmon2", "mediumpurple4", "orange", "maroon", "yellow3", "brown4", "yellow4", "sienna4", "chocolate", "gray19")
# FTE Project
fte_merged_all_filt <- readRDS('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/ReProcessed/Analysis_RDS/integrated_bySample_filtered.rds')


plotCellTypeProps(clusters=fte_merged_all_filt$CellType_manual, sample=fte_merged_all_filt$Mutation) + scale_fill_manual(values = kelly.colours[1:18])


### Propeller test 
output.asin <- propeller(clusters=fte_merged_all_filt$CellType_manual, sample=fte_merged_all_filt$orig.ident, 
                         group=fte_merged_all_filt$Mutation, transform="logit")
output.asin


props.fte <- getTransformedProps(clusters=fte_merged_all_filt$CellType_manual, sample=fte_merged_all_filt$orig.ident,
                                   transform="logit")

### Match sample id with Mutation status
# fte_meta <- fread('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/Tables/Metadata_AllSamples.tsv')
# fte_meta$ID <- gsub(' ', '', fte_meta$ID)
# fte_meta$ID <- paste0('sample_',fte_meta$ID)
# 


plot.prop <- as.data.frame(props.fte$Proportions)
plot.prop$Mutation <- as.factor(plot.prop$sample)
levels(plot.prop$Mutation) <- c('BRCA2', 'WT', 'BRCA2', 'MLH1', 'MLH1', 
                                'PALB2', 
                                'BRCA1', 
                                'BRCA1', 
                                'WT',
                               'BRCA2',
                               'BRCA2', 
                               'BRCA2',
                               'BRCA2', 'WT',
                               'WT', 'PALB2',
                               'Tumor/Tumor+BRCA1', 'PALB2', 'Tumor/Tumor+BRCA1', 
                               'RAD51C', 'Tumor/Tumor+BRCA1', 'BRCA1',
                               'PALB2', 
                               'BRCA2', 
                               'BRCA2', 
                               'BRCA1',
                               'BRCA1',
                               'RAD51C', 'RAD51C', "WT")
plot.prop$Mutation <- as.character(plot.prop$Mutation)


tag_facet <- function(p, open = "(", close = ")", tag_pool = letters, x = -Inf, y = Inf, 
                      hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
  
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
  p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE) 
}
p <- ggplot(plot.prop, aes(Mutation, Freq, color = Mutation)) + geom_point(size=2.5) + facet_wrap(~clusters) +
  xlab('Mutation') + ylab('Proportions') + theme(text = element_text(size = 17)) + coord_flip() 

my_tag <- c('p= 0.247', 'p= 0.0737', 'p= 0.666', 'p= 0.010',
            'p= 0.666', 'p= 0.284', 'p= 0.466', 'p= 0.055', 'p= 0.261',
            'p= 0.010', 'p= 0.910', 'p= 0.060', 'p= 0.102', 'p= 0.910',
            'p= 0.038', 'p= 0.248', 'p= 0.613', 'p= 0.985')
my_tag <- paste0('p= ', prettyNum(output.asin[levels(as.factor(plot.prop$clusters)), ]$FDR, digits=2) )

tag_facet(p, 
          x = -Inf, y = -Inf, 
          vjust = -2, hjust = -2,
          open = "", close = "",
          fontface = 1,
          size = 4,
         # family = "serif",
          tag_pool = my_tag)



### Monocle 
library(SeuratWrappers)
library(monocle3)
cds <- as.cell_data_set(fte_merged_all_filt)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "CellType")

cds <- cluster_cells(cds, resolution=1e-3)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "CellType",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds <- order_cells(cds, reduction_method = 'UMAP')
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

                       




