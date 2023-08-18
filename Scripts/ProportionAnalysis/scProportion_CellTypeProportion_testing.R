#Single Cell Proportion Test
library(scProportionTest)


prop_test <- sc_utils(ft1)

prop_test <- permutation_test(
  prop_test, cluster_identity = "CellType_manual",
  sample_1 = "WT", sample_2 = "BRCA1",
  sample_identity = "Mutation"
)
permutation_plot(prop_test)


## For all mutations vs WT

mutations <- levels(as.factor(ft1$Mutation)) %>% as.character
# remove WT since we dont want to compare WT with WT lol
mutations <- mutations[-7]

for(i in mutations) {
  prop_test <- permutation_test(
    prop_test, cluster_identity = "CellType_manual",
    sample_1 = "WT", sample_2 = i,
    sample_identity = "Mutation"
  )
  p = permutation_plot(prop_test) + ggtitle(paste0('WT vs ', i)) + theme(text = element_text(size = 16))  
  ggsave(paste0('/mnt/plummergrp/Felipe/FTE_scRNA_Miami/plots/PermutationPlot_CellTypeProportions_WT_vs_', i, '.png'), plot = p, width = 10, height = 10)
}
