##### Dr George scRNA Project -- FTE tissue; ######
### RE PROCESSED DATA; cellranger 7.1; intron mode on. ; July 2023 - Felipe Segato Dezem ###





## Pre-processing - run SoupX and Standard QC filtering using Seurat
library(SoupX)
library(Seurat)



master_dir <- '/media/ResearchHome/plummgrp/projects/SpatialAssays_Jiyang/common/bkp_felipe/cellranger_outs/'




tmp.dir <- dir(master_dir)
seurat_qc_table <- data.frame(SampleID = gsub("_Cell.*", '', tmp.dir), nCells = NA, nCells_filt = NA,
                              nGenes = NA, nGenes_filt = NA)
for(i in 1:length(tmp.dir)) {
  
  ## Setting the correct folder name according to the output structure of different versions of CellRanger. zzz chato pra caralho zzz
  folders_in <- dir(paste0(master_dir, tmp.dir[i],'/outs'))[grep("filtered_[feature|gene]", dir(paste0(master_dir, tmp.dir[i],'/outs')))]
  ##remove the .h5 file and leave only the folder name 
  w <- which(folders_in %in% c("filtered_gene_bc_matrices" , "filtered_feature_bc_matrices", "filtered_gene_bc_matrix", 
                               "filtered_feature_bc_matrix"))
  folders_in <- folders_in[w]
  
  ##Condition because depending on which version theres a subfolder 
  if(folders_in == "filtered_feature_bc_matrix") {
    sc = Seurat::Read10X(paste0(master_dir, tmp.dir[i], '/outs/', folders_in, "/"))
  } else {
    sc = Seurat::Read10X(paste0(master_dir, tmp.dir[i], '/outs/', folders_in, "/GRCh38/"))
    
  }
  
  sc = CreateSeuratObject(sc, min.cells = 3, min.features = 200, 
                          project = gsub("_Cell.*", '', tmp.dir[i]))
  # update qc table with # of cells pre-filt
  seurat_qc_table[i, ]$nCells <- dim(sc)[2]; seurat_qc_table[i, ]$nGenes <- dim(sc)[1];
  #
  
  
  # Filter poor quality cells
  sc[['percent.mt']] = PercentageFeatureSet(sc, pattern = '^MT-')
  sc[['percent.ribo']] = PercentageFeatureSet(sc, pattern ="^RP[SL]")
  #VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  sc = subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25 & nCount_RNA > 200)
  
  # update qc table with # of cells aftter filtering
  seurat_qc_table[i, ]$nCells_filt <- dim(sc)[2]; seurat_qc_table[i, ]$nGenes_filt <- dim(sc)[1];
  #
  
  # Cluster 
  # sc = NormalizeData(sc, normalization.method = 'LogNormalize', scale.factor = 10000, verbose = FALSE)
  # sc = FindVariableFeatures(sc, selection.method = 'vst', nfeatures = 2000, verbose = FALSE)
  # sc = ScaleData(sc, verbose = FALSE, features = rownames(sc))
  # sc = RunPCA(sc, verbose = FALSE)
  # sc = RunUMAP(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
  # sc = FindNeighbors(sc, reduction = 'pca', dims = 1:20, verbose = FALSE)
  # sc = FindClusters(sc, resolution = 0.5, verbose = FALSE)
  
  
  
  # ### SoupX
  # toc = Seurat::Read10X(file.path(paste0(master_dir, tmp.dir[i], '/outs/'), "filtered_gene_bc_matrices", "GRCh38"))
  # tod = Seurat::Read10X(file.path(paste0(master_dir, tmp.dir[i], '/outs/'), "raw_gene_bc_matrices", "GRCh38"))
  # soup = SoupChannel(tod, toc)
  # soup = setClusters(soup, sc@active.ident)
  # soup = autoEstCont(soup)
  # out = adjustCounts(soup)
  # srat = CreateSeuratObject(out)
  # # p1 <- VlnPlot(srat, features = c('nCount_RNA', 'nFeature_RNA'))
  # # p2 <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
  # # p2 + p1
  
  # Save data as RDS file 
  saveRDS(sc, file = paste0('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/ReProcessed/', 
                            gsub("_Cell.*", '', tmp.dir[i]), '_filtered.rds'))
  
  print(paste0('Done processing sample ', tmp.dir[i]))
  
  
}

seurat_qc_table$Percent_filtered <- prettyNum(((seurat_qc_table$nCells-seurat_qc_table$nCells_filt) / seurat_qc_table$nCells) * 100, digits =2) %>%
  as.numeric()
mean(seurat_qc_table$Percent_filtered)

write.table(seurat_qc_table, row.names=F, quote=F, file = '/mnt/plummergrp/Felipe/FTE_scRNA_Miami/Tables/ReProcessed_QC_Table.tsv',
            sep = "\t")



