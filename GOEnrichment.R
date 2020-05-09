## R 3.6.3
## -*- coding: utf-8 -*-
## Gene Enrichment Project.

# Packages preperation.
install.packages('devtools', repos = 'http://mirrors.tuna.tsinghua.edu.cn/CRAN')
install.packages('BiocManager', repos = 'http://mirrors.tuna.tsinghua.edu.cn/CRAN')
install.packages('igraph')
library(devtools)
BiocManager::install(version = '3.10')
Needed <- c('bit', 'formatR', 'hms', 'triebeard', 'tweenr', 'polyclip',
           'RcppEigen', 'RcppArmadillo', 'zlibbioc', 'bit64', 'blob', 
           'plogr', 'lambda.r', 'futile.options', 'progress', 'urltools',
           'gridGraphics', 'ggforce', 'ggrepel', 'viridis', 'tidygraph',
           'graphlayouts', 'bitops', 'XVector', 'IRanges', 'RSQLite',
           'futile.logger', 'snow', 'data.table', 'gridExtra', 'fastmatch',
           'cowplot', 'europepmc', 'ggplotify', 'ggraph', 'ggridges','igraph', 
           'dplyr', 'tidyselect', 'RCurl', 'Biostrings', 'AnnotationDbi', 
           'BiocParallel', 'DO.db', 'fgsea', 'GOSemSim','qvalue', 'S4Vectors', 
           'BiocGenerics', 'graph', 'Biobase','GO.db', 'SparseM', 'matrixStats', 
           'DBI', 'enrichplot', 'rvcheck', 'tidyr', 'org.Hs.eg.db', 'KEGGgraph', 
           'XML', 'Rgraphviz', 'png', 'KEGGREST', 'GOplot')
install.packages(Needed, repos = 'http://mirrors.tuna.tsinghua.edu.cn/CRAN')
BiocManager::install(c('DOSE', 'topGO', 'clusterProfiler', 'pathview'))

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(igraph)
library(ggplot2)
library(stringr)
library(GOplot)

# ---------------------------------------------------------------------------------
# To load target gene data.
# Method 1:
MyGeneSet <-c('TP53', 'TNF', 'EGFR', 'ESR1', 'BRCA1', 'APC')
MyGeneIDSet = bitr(MyGeneSet,
                   fromType = 'SYMBOL',
                   toType = 'ENTREZID',
                   OrgDb = 'org.Hs.eg.db')

# Method 2:
setwd('E:/Learning/2019_09-2020_08/Spring/Language/Task_GO')
MyGeneSet2 <- read.table('GeneNames.txt', header = T)
MyGeneSet2$SYMBOL <- as.character(MyGeneSet2$SYMBOL)
MyGeneIDSet2 = bitr(MyGeneSet2,
                   fromType = 'SYMBOL',
                   toType = 'ENTREZID',
                   OrgDb = 'org.Hs.eg.db')
# To load the whole gene list.
data(geneList, package = 'DOSE')


# GO enrichment: all.
ego_ALL <- enrichGO(gene = MyGeneSet2$ENTREZID,
                    universe = names(geneList),
                    OrgDb = org.Hs.eg.db,
                    ont = 'ALL',
                    pAdjustMethod = 'BH',
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = T)
write.csv(summary(as.data.frame(ego_ALL)), 'GO-Enrich-ALL.csv', row.names = F)
# GO enrichment: Molecular Function.
ego_MF <- enrichGO(gene = MyGeneSet2$ENTREZID,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   ont = 'MF',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = T)
# GO enrichment: Cell Component.
ego_CC <- enrichGO(gene = MyGeneSet2$ENTREZID,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   ont = 'CC',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = T)
# GO enrichment: Biological Process.
ego_BP <- enrichGO(gene = MyGeneSet2$ENTREZID,
                   universe = names(geneList),
                   OrgDb = org.Hs.eg.db,
                   ont = 'BP',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1,
                   readable = T)


# Enrichment analysis and visualization.
ego_MF_df <- as.data.frame(ego_MF)
View(ego_MF_df) # 83
ego_CC_df <- as.data.frame(ego_CC)
View(ego_CC_df) # 87
ego_BP_df <- as.data.frame(ego_BP)
View(ego_BP_df) #1495
#ego_MF_df <- ego_MF_df[order(ego_MF_df[ , 'Count'], decreasing = T), ]
#ego_CC_df <- ego_CC_df[order(ego_CC_df[ , 'Count'], decreasing = T), ]
#ego_BP_df <- ego_BP_df[order(ego_BP_df[ , 'Count'], decreasing = T), ]

# An overall display of GO enrichment result.
show_num <- c(25, 25, 50)
ego_result_MF <- ego_MF_df[1: show_num[1], ]
ego_result_CC <- ego_CC_df[1: show_num[2], ]
ego_result_BP <- ego_BP_df[1: show_num[3], ]
go_result_df <- data.frame(ID = c(ego_result_BP$ID, 
                                  ego_result_CC$ID, 
                                  ego_result_MF$ID),
                            Description = c(ego_result_BP$Description, 
                                            ego_result_CC$Description, 
                                            ego_result_MF$Description),
                            Count = c(ego_result_BP$Count, 
                                      ego_result_CC$Count, 
                                      ego_result_MF$Count),
                            GO_type = factor(c(rep('BP', length(ego_result_BP$ID)), 
                                               rep('CC', length(ego_result_CC$ID)), 
                                               rep('MF', length(ego_result_MF$ID)))))
go_result_df$Serial_num <- factor(rev(1:nrow(go_result_df)))
cols <- c('#8DA1CB', '#FD8D62', '#66C3A5')
x_label <- levels(go_result_df$Description)
names(x_label) <- rev(1: nrow(go_result_df))
ggplot(data = go_result_df, aes(x = Serial_num, y = Count, fill = GO_type)) +
  geom_bar(stat = "identity", width = 0.8) + coord_flip() + 
  scale_fill_manual(values = cols) + theme_bw() + 
  scale_x_discrete(labels = x_label) +
  xlab('GO term') +
  theme(axis.text = element_text(face = "bold", color="gray40", size = 7)) +
  labs(title = "General Display of Enriched GO Terms")

# Specific statistical visualization: GO - ALL.
par(cex = 0.45)
barplot(ego_ALL, showCategory = 20, title = 'GO Enrichment Barplot - ALL')
dotplot(ego_ALL, showCategory = 50, title = 'GO Enrichment Dotplot - ALL')
cnetplot(ego_ALL, foldChange = geneList, title = 'GO Enrichment Netplot - ALL')
# Specific statistical visualization: GO - MF.
barplot(ego_MF, showCategory = 20, title = 'GO Enrichment Barplot - MF')
dotplot(ego_MF, showCategory = 50, title = 'GO Enrichment Dotplot - MF')
cnetplot(ego_MF, foldChange = geneList, title = 'GO Enrichment Netplot - MF')
plotGOgraph(ego_MF)
# Specific statistical visualization: GO - CC.
barplot(ego_CC, showCategory = 20, title = 'GO Enrichment Barplot - CC')
dotplot(ego_CC, showCategory = 50, title = 'GO Enrichment Dotplot - CC')
cnetplot(ego_CC, foldChange = geneList, title = 'GO Enrichment Netplot - CC')
plotGOgraph(ego_CC)
# Specific statistical visualization: GO - BP.
barplot(ego_BP, showCategory = 20, title = 'GO Enrichment Barplot - BP')
dotplot(ego_BP, showCategory = 50, title = 'GO Enrichment Dotplot - BP')
cnetplot(ego_BP, foldChange = geneList, title = 'GO Enrichment Netplot - BP')
plotGOgraph(ego_BP)


# To rebuild the file format that can be recognized by GOplot.
GO <- ego_ALL[1: 20, c(1, 2, 3, 9, 7)]
GO$geneID <- str_replace_all(GO$geneID, "/", ",")
names(GO) <- c("Category","ID", "term", "Genes", "adj_pval")
gene <- data.frame(ID = MyGeneIDSet$SYMBOL, logFC = rnorm(length(MyGeneIDSet$ENTREZID), mean=0, sd=2))
# To use GOplot for visualization.
circ <- circle_dat(GO,gene)
GOBar(circ)
GOCircle(circ)
chord <- chord_dat(data = circ, genes = gene, process = GO$term)
GOChord(chord,space = 0.02, gene.order = "logFC", gene.space = 0.25, gene.size = 5)
GOHeat(chord[, 1: 10], nlfc = 0)#nlfc=0 no logFC


# KEGG enrichment analysis.
kegg <- enrichKEGG(gene = MyGeneIDSet$ENTREZID, universe = names(geneList), 
                   organism ="hsa", pvalueCutoff = 0.01, pAdjustMethod = 'BH',
                   minGSSize = 3, maxGSSize = 500, qvalueCutoff = 0.01,
                   use_internal_data = FALSE)

barplot(kegg, showCategory = 20, title = 'KEGG Enrichment Bar')
dotplot(kegg, showCategory = 50, title = 'KEGG Enrichment Dot')
cnetplot(kegg, foldChange = geneList, title = 'KEGG Enrichment Net')

