#### SCRIPT JUST FOR PLOTTING FIGURES FOR THE FT PAPER #######



## UMAP no tumors 
pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/MainUMAP_CellType_NoTumors.pdf', width =20, height = 10)
DimPlot(ft1, reduction = 'umap.harmony', cells = colnames(ft1)[!ft1$Mutation %in% c( 'Tumor/Tumor+BRCA1')], 
        group.by = c('CellType_manual'), label = T) & NoAxes()
dev.off()


pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/BarPlot_nCells_bySample_byMutation.pdf', width =20, height = 10)
dittoBarPlot(
  object = ft1[, !ft1$Mutation %in% c( 'Tumor/Tumor+BRCA1') ],
  var = "orig.ident",
  group.by = "Mutation", scale ='count') 
dev.off()


## UMAP no tumors 
pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/UMAP_Integration_NoTumors.pdf', width =30, height = 10)
DimPlot(ft1, reduction = 'umap.harmony', cells = colnames(ft1)[!ft1$Mutation %in% c( 'Tumor/Tumor+BRCA1')], 
        group.by = c('orig.ident', 'Mutation', 'harmony_clusters_res0.5')) & NoAxes()
dev.off()


## UMAP Clusters split by Mutation
pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/UMAP_CellType_byMutation_NoTumors.pdf', width =30, height = 7)
ft1.sub <- ft1[, !ft1$Mutation %in% c( 'Tumor/Tumor+BRCA1')]
DimPlot(ft1.sub, reduction = 'umap.harmony', 
        group.by = c('CellType_manual'), split.by = 'Mutation') & NoAxes()
dev.off()


## Barplot proportion celltype per sample / mutation


pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/Barplot_Proportion_Celltype_byMutation.pdf', width =20, height = 10)
dittoBarPlot(
  object = ft1.sub,
  var = "CellType_manual",
  group.by = "Mutation") 
dev.off()



pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/Barplot_Proportion_Celltype_bySample.pdf', width =20, height = 10)
dittoBarPlot(
  object = ft1.sub,
  var = "CellType_manual",
  group.by = "orig.ident") 
dev.off()


## Dotplot markers 
markers2 <- c('EPCAM', 'FOXJ1', 'CAPS', 'KRT7', 'PAX8', 'OVGP1', 'DCN', 'COL1A1', 'CD34', 'PDGFRA',
              'POSTN', 'NR2F2', 'DES', 'ACTA2', 'MYH11', 'TAGLN', 'PDGFRB', 'MCAM', 'CSPG4', 'KDR', 'PECAM1',
              'VWF', 'PROX1', 'PDPN', 'SDC1', 'JCHAIN', 'MS4A1', 'KLRC1', 'RUNX3', 'CD3E', 'PTPRC', 'KIT',
              'MS4A2', 'TPSB2', 'TPSAB1', 'FOLR2', 'CD68', 'CD163', 'ITGAX')
pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/DotPlot_CellType_markers.pdf', width =20, height = 10)
DotPlot(ft1.sub, features = markers2, group.by = 'CellType_manual') & RotatedAxis() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()




## FeaturePlot of markers
pdf('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/FigurePlots/FeaturePlot_CellType_markers.pdf', width =50, height = 50)
FeaturePlot(ft1.sub, features = markers2, reduction = 'umap.harmony') & NoAxes() & scale_color_viridis()
dev.off()


# dotplot for the general differences in Mutation by celltype was plotted on the propeller script.












