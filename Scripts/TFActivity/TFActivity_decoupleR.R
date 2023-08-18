## devcoupleR scRNA Ft data

library(decoupleR)
library(OmnipathR)

# Extra libraries
library(dplyr)
library(pheatmap)

net <- get_collectri(organism='human', split_complexes=FALSE)


# Extract the normalized log-transformed counts
mat <- as.matrix(ft1@assays$RNA$data)

# Run ulm
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)

# Extract ulm and store it in tfsulm in pbmc
ft1[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)

# Change assay
DefaultAssay(object = ft1) <- "tfsulm"

# Scale the data
ft1 <- ScaleData(ft1)
ft1@assays$tfsulm@data <- ft1@assays$tfsulm@scale.data

p1 <- DimPlot(ft1, reduction = "umap.harmony", label = TRUE, pt.size = 0.5) + 
  NoLegend() + ggtitle('Cell types')
p2 <- (FeaturePlot(ft1, features = c("PAX8"), reduction = 'umap.harmony') & 
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')) +
  ggtitle('PAX8 activity')
DefaultAssay(ft1) <- "RNA"
p3 <- FeaturePlot(ft1, features = c("PAX8"), reduction = 'umap.harmony') + ggtitle('PAX8 expression') + scale_color_viridis()
DefaultAssay(object = ft1) <- "tfsulm"
p1 | p2 | p3



## get top 50 TFs per Mutation type 
n_tfs <- 50
mutation <- levels(as.factor(ft1$Mutation)) %>% as.character
df_list <- list()
tfs_list <- list()
top_acts_mat_list <- list()

for(i in mutation) {
  # Extract activities from object as a long dataframe
  df_list[[i]] <- t(as.matrix(ft1[, ft1$Mutation %in% i]@assays$tfsulm@scale.data)) %>%
    as.data.frame() %>%
    mutate(cluster = Idents(ft1[, ft1$Mutation %in% i])) %>%
    pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
    group_by(cluster, source) %>%
    summarise(mean = mean(score))
  
  # Get top tfs with more variable means across clusters
  tfs_list[[i]] <- df_list[[i]] %>%
    group_by(source) %>%
    summarise(std = sd(mean)) %>%
    arrange(-abs(std)) %>%
    head(n_tfs) %>%
    pull(source)
  
  # Subset long data frame to top tfs and transform to wide matrix
  top_acts_mat_list[[i]] <- df_list[[i]] %>%
    filter(source %in% tfs_list[[i]]) %>%
    pivot_wider(id_cols = 'cluster', names_from = 'source',
                values_from = 'mean') %>%
    column_to_rownames('cluster') %>%
    as.matrix()
  
}



# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-4, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 4, length.out=floor(palette_length/2)))

# Plot
for(i in mutation) {
  p = pheatmap(top_acts_mat_list[[i]], border_color = NA, color=my_color, breaks = my_breaks, main = i) 
  ggsave(paste0('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/plots/HM_TFActivity_', i, '.png'),plot =p, width = 10, height = 8)
  
}


pheatmap(top_acts_mat, border_color = NA, color=my_color, breaks = my_breaks, main = 'WT') 


### DorotheA Regulons
library(dorothea)

## We read Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## We obtain the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## We compute Viper Scores 
DefaultAssay(ft1) <- 'RNA'
ft1 <- run_viper(ft1, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 10, 
                                 verbose = T))

## We compute the Nearest Neighbours to perform cluster
DefaultAssay(object = ft1) <- "dorothea"
ft1 <- ScaleData(ft1)
ft1 <- RunPCA(ft1, features = rownames(ft1), verbose = FALSE)
ft1 <- FindNeighbors(ft1, dims = 1:10, verbose = FALSE)
ft1 <- FindClusters(ft1, resolution = 0.5, verbose = FALSE, cluster.name = 'regulon_clusters')
ft1 <- RunUMAP(ft1, dims = 1:10, umap.method = "uwot", metric = "cosine", reduction.key = 'umap.regulon')
#

## We transform Viper scores, scaled by seurat, into a data frame to better 
## handling the results
viper_scores_df <- GetAssayData(ft1, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame() %>%
  t()

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = names(Idents(ft1)), 
                            cell_type = as.character(ft1$CellType_manual),
                            stringsAsFactors = FALSE)
CellsClusters$cell <- gsub('\\-', '.', CellsClusters$cell)

## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))


# We select the 20 most variable TFs. (20*9 populations = 180)
highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(500, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "DoRothEA (ABC)", angle_col = 45,
                       treeheight_col = 0,  border_color = NA) 

