### Pathway ORA on cluster markers scRNA ft 

library(clusterProfiler)
library(tidyr)
library(org.Hs.eg.db)
library(ggplot2)

degs <- read.table('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/Tables/CellType_top50Markers_clusters_res0.5.csv', sep = '\t' ,header = T)

cell_types = levels(as.factor(degs$group))

ego_list <- list()
for(i in cell_types) {
    
  ego <- enrichGO(gene         =  degs[degs$group %in% i, ]$feature,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  
  ego_list[[i]] <- ego

  
  
}


for(i in 1:length(ego_list)) {
  p = dotplot(ego_list[[i]], showCategory=20) & scale_color_viridis() & ggtitle(names(ego_list)[i])
  ggsave(paste0('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/plots/dotplot_ORA_celltype_', i, '.png'), plot = p, width = 10, height = 16)
}
