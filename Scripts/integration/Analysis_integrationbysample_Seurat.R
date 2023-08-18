
# We reviewed the clusters as a group but had difficulty drawing out the biological story with just the raw gene lists for each cluster. Could you provide us with a few other things to assist with this?
#   
#   Gene ontology/ GSEA pathway analysis for each cluster
# Differential expression analysis comparing each of the ciliated epithelia clusters. Also for the secretory epithelia clusters 
# UMAP with just the tumors --> highly expressed genes (that we can then search for in the mutation carriers)
# Some other UMAPs that are subsampled - one for all epithelial cell types, one for fibroblasts, one for all immune cells. And then the relevant differential expression and pathway analysis
# Is the EMT-like a ciliated or secretory epithelial cluster?
#   Also of note, all the tumors and MLH1 samples should be excluded from this initial plot
# 
# After the getting this additional data, we would like to do a few additional analyses for the paper: ligand-receptor interaction, Cell state scoring, immune enrichment scoring, T cell exhaustion. These can be done using a few different R packages. Do you have time/capacity to do this? If not we can also pursue this with our bioinformaticians. 
# 



## Load saved object

ft1 = readRDS('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/ReProcessed/Analysis_RDS/integrated_bySample_filtered.rds')
ft1 <- JoinLayers(ft1)

## Update metadata with Mutation information
ft1$Mutation <- NA
ft1@meta.data[ft1$orig.ident %in% c('sample_199', 'sample_EOC70', 'sample_5010', 'sample_5137', 'sample_EOC75_merged', 'sample_FTE1222',
                                                          'sample_FTE1322', 'sample_FTE1722') , ]$Mutation <- 'BRCA1'

ft1@meta.data[ft1$orig.ident %in% c('sample_3313', 'sample_FTE1122', 'sample_5334', 'sample_5334_11713', 'sample_5349', 'sample_FTE0721', 'sample_4586',
                                                          'sample_5273', 'sample_5334') , ]$Mutation <- 'BRCA2'

ft1@meta.data[ft1$orig.ident %in% c('sample_EOC72_merged', 'sample_FTE1422', 'sample_FTE1522') , ]$Mutation <- 'RAD51C'

ft1@meta.data[ft1$orig.ident %in% c('sample_4630', 'sample_4794') , ]$Mutation <- 'MLH1'

ft1@meta.data[ft1$orig.ident %in% c('sample_EOC56_merged', 'sample_4931', 'sample_EOC27', 'sample_FTE0521', 'sample_FTE0821') , ]$Mutation <- 'PALB2'

#ft1@meta.data[ft1$orig.ident %in% c('sample_EOC53', 'sample_EOC55') , ]$Mutation <- 'Tumor'

ft1@meta.data[ft1$orig.ident %in% c('sample_EOC53', 'sample_EOC55','sample_EOC71', 'sample_EOC74') , ]$Mutation <- 'Tumor/Tumor+BRCA1'

ft1@meta.data[ft1$orig.ident %in% c('sample_192', 'sample_5220', 'sample_BSSR2676', 'sample_FTE1822', 'sample_4354', 'sample_BSSR2896') , ]$Mutation <- 'WT'
###





## cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

ft1 <- CellCycleScoring(ft1, s.features = s.genes, g2m.features = g2m.genes)

DotPlot(ft1, features = markers2, group.by = 'customclassif')  & RotatedAxis() &  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))



## Semi supervised cell annotation 
# lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
# # load gene set preparation function
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# # load cell type annotation function
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
# 
# db_ = '/mnt/plummergrp/maycon/Miami_Ovarian_scrnaseq/high_intersect_93_markers_low_intersect_18_markers_no_orf_gene.xlsx'
# db_ = '/mnt/plummergrp/Felipe/FTE_scRNA_Miami/Tables/db_celltype_david_updated.xlsx'
# tissue = "FT"
# 
# 
# # prepare gene sets
# library(stringr)
# 
# cell_markers = openxlsx::read.xlsx(db_)
# gs_list <- strsplit(cell_markers$geneSymbolmore1, ',')
# names(gs_list) <- cell_markers$cellName
# gs_list <- lapply(gs_list, function(x) {str_trim(x, "left") })
# 
# 
# # get cell-type by cell matrix
# es.max = sctype_score(scRNAseqData = ft1[["RNA"]]$data, scaled = F, 
#                       gs = gs_list, gs2 = NULL)
# 
# # merge by cluster
# cL_resutls = do.call("rbind", lapply(unique(ft1@meta.data$harmony_clusters_res0.5), function(cl){
#   es.max.cl = sort(rowSums(es.max[ ,rownames(ft1@meta.data[ft1@meta.data$harmony_clusters_res0.5==cl, ])]), decreasing = !0)
#   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ft1@meta.data$harmony_clusters_res0.5==cl)), 20)
# }))
# sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 2, wt = scores)  
# 
# # set low-confident (low ScType score) clusters to "unknown"
# sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unknown"
# print(sctype_scores[,1:3])
# 
# ft1@meta.data$customclassif = ""
# for(j in unique(sctype_scores$cluster)){
#   cl_type = sctype_scores[sctype_scores$cluster==j,]; 
#   ft1@meta.data$customclassif[ft1@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
# }
# DimPlot(ft1, group.by = c('customclassif','orig.ident', 'Mutation'), label=T, reduction = 'umap.harmony')
#### Didn't work very wel... - trying to do it manually



## Markers fromn paper FT development
library(viridis)
markers2 <- c('EPCAM', 'FOXJ1', 'CAPS', 'KRT7', 'PAX8', 'OVGP1', 'DCN', 'COL1A1', 'CD34', 'PDGFRA',
              'POSTN', 'NR2F2', 'DES', 'ACTA2', 'MYH11', 'TAGLN', 'PDGFRB', 'MCAM', 'CSPG4', 'KDR', 'PECAM1',
              'VWF', 'PROX1', 'PDPN', 'SDC1', 'JCHAIN', 'MS4A1', 'KLRC1', 'RUNX3', 'CD3E', 'PTPRC', 'KIT',
              'MS4A2', 'TPSB2', 'TPSAB1', 'FOLR2', 'CD68', 'CD163', 'ITGAX')
DotPlot(ft1, features = markers2, group.by = 'harmony_clusters_res0.5') & RotatedAxis() & scale_color_viridis()

## Annotate based on these markers - manual
ft1$CellType_manual <- as.factor(ft1$harmony_clusters_res0.5) 
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('9')] <- 'CD4 TCell'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('8')] <- 'CD8 TCell'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('15')] <- 'NK Cell'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('10', '7')] <- 'Ciliated Epithelial'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('5')] <- 'Non-ciliated Secretory Epithelial 1'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('0','24')] <- 'Non-ciliated Secretory Epithelial 2'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('2', '23')] <- 'Non-ciliated Secretory Epithelial 3'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('3')] <- 'Non-ciliated Secretory Epithelial 4'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('22')] <- 'Mast'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('19')] <- 'BCell'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('13')] <- 'Blood Endothelial'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('12')] <- 'Pericyte'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('4')] <- 'Smooth Muscle'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('23')] <- 'Myofibroblast'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('1', '27')] <- 'Fibroblast'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('20')] <- 'Fibroblast, DCN-'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('14')] <- 'EMT-like'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('6')] <- 'Macrophage 1'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('11')] <- 'Macrophage 2'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('17')] <- 'Epithelial, EPCAM+ OVGP1+ PAX8-'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('28')] <- 'Fibroblast, S100A11+'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('25')] <- 'Epithelial, EPCAM+ CAPS+ OVGP1+ PAX8-'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('16')] <- 'Myeloid'
levels(ft1$CellType_manual)[levels(ft1$CellType_manual) %in% c('18')] <- 'Endothelial'

  library(viridis); library(RColorBrewer)
DimPlot(ft1, group.by = 'CellType_manual', label = T, reduction = 'umap.harmony')
DotPlot(ft1, features = markers2, group.by = 'CellType_manual') & RotatedAxis() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


####
sub1 <- subset(ft1, subset = Mutation %in% c('WT', 'Tumor/Tumor+BRCA1'))




### Comparing epithelial clusters cilliated / non-ciliated ##

# Idents(ft1) <- ft1$CellType_manual
# deg_secepi1 <- FindMarkers(ft1, ident.1 = 'Non-ciliated Secretory Epithelial 1', ident.2 = 'Non-ciliated Secretory Epithelial 2', logfc.threshold = 1)
# deg_secepi2 <- FindMarkers(ft1, ident.1 = 'Non-ciliated Secretory Epithelial 1', ident.2 = 'Non-ciliated Secretory Epithelial 3', logfc.threshold = 1)
# 

## Subsetting the object with the cell types of interest and then performing the one-vs-all rank-sum wilcoxon test
library(presto)
# Non ciliated
degs_noncil <- presto::wilcoxauc(subset(ft1, CellType_manual %in% c(paste0('Non-ciliated Secretory Epithelial ', seq(1:4)))), 'CellType_manual')
top_degs_noncil = degs_noncil %>% group_by(group) %>% top_n(n = 20, wt = logFC)

DotPlot(subset(ft1, CellType_manual %in% c(paste0('Non-ciliated Secretory Epithelial ', seq(1:4)))),
        features = unique(top_degs_noncil$feature), group.by = 'CellType_manual') & RotatedAxis() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

# All epithelial cells
degs_epit <- presto::wilcoxauc(subset(ft1, CellType_manual %in% c(paste0('Non-ciliated Secretory Epithelial ', seq(1:4)), 'Ciliated Epithelial', 'Epithelial, EPCAM+ CAPS+ OVGP1+ PAX8-', 'Epithelial, EPCAM+ OVGP1+ PAX8-')), 
                               'CellType_manual')

top_epi = degs_epit %>% group_by(group) %>% top_n(n = 10, wt = logFC)

DotPlot(subset(ft1, CellType_manual %in% c(paste0('Non-ciliated Secretory Epithelial ', seq(1:4)), 'Ciliated Epithelial', 'Epithelial, EPCAM+ CAPS+ OVGP1+ PAX8-', 'Epithelial, EPCAM+ OVGP1+ PAX8-')),
        features = unique(top_epi$feature), group.by = 'CellType_manual') & RotatedAxis() & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


## UMAP with only tumors

DimPlot(ft1, reduction = 'umap.harmony', cells = colnames(ft1)[ft1$Mutation %in% 'Tumor/Tumor+BRCA1'], group.by = c('CellType_manual','Mutation', 'orig.ident')) & NoAxes()


## UMAP with only epithelial cells

DimPlot(ft1, reduction = 'umap.harmony', cells = colnames(ft1)[ft1$CellType_manual %in% c(paste0('Non-ciliated Secretory Epithelial ', seq(1:4)), 'Ciliated Epithelial', 'Epithelial, EPCAM+ CAPS+ OVGP1+ PAX8-', 'Epithelial, EPCAM+ OVGP1+ PAX8-')], 
        group.by = c('CellType_manual','Mutation', 'orig.ident')) & NoAxes()

## UMAP fibroblasts
DimPlot(ft1, reduction = 'umap.harmony', cells = colnames(ft1)[ft1$CellType_manual %in% c('Fibroblast, S100A11+', 'Fibroblast', 'Fibroblast, DCN-')], 
        group.by = c('CellType_manual','Mutation', 'orig.ident')) & NoAxes()

## Immune cells
DimPlot(ft1, reduction = 'umap.harmony', cells = colnames(ft1)[ft1$CellType_manual %in% c('Macrophage 2', 'Mast', 'CD8 TCell', 'Pericytes', 'NK Cell', 'CD4 TCell', 'BCell', 'Macrophage 1')], 
        group.by = c('CellType_manual','Mutation', 'orig.ident')) & NoAxes()

















