library(data.table)
library(openxlsx)
library(mygene)
library(dplyr)
library(tidygraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(grid)
library(igraph)

source('/share/home/sunly/Rscript/cbmap_pro/cbmap_586/2_sample_process.R')
adnc_lim <- read.xlsx('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/table/dea_result/cbmap_adnc_lim_re.xlsx')
ppi <- as.data.frame(fread('/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/processed/processed_data/wmds/cbmap_586/physical_9/ppi_select.txt'))
dep_list <- read.xlsx('/data/projects/China_Brain_MultiOmics/humanBrain_protein/WMDSnet/input/cbmap_586_stringPhysical9_LMH_regCOV_WMDSnet_revise_geneWeight_list.xlsx')
enrich_kegg <- as.data.frame(fread('/data/projects/China_Brain_MultiOmics/humanBrain_protein/WMDSnet/result/cbmap_586/WMDS.net_KEGG_20250930.txt')[, -1])
ad_gene_literature <- read.csv('alzheimer_gene_literature_final.csv')

ad_gene_literature_select <- ad_gene_literature[-which(ad_gene_literature$type=='unknown'),]
ad_gene_num <- ad_gene_literature_select %>%
  count(gene_name, name = "num")

ad_gene <- enrich_kegg %>%
  filter(Description == "Alzheimer disease") %>%
  pull(geneID_symbols) %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

ppi_select <- ppi %>%
  filter(
    gene_A %in% dep_list$Gene,
    gene_B %in% dep_list$Gene
  )

all_nodes <- unique(c(ppi_select$gene_A, ppi_select$gene_B))
node <- data.frame(node = all_nodes)
node$Weight <- dep_list$Weight[match(node$node, dep_list$Gene)]
node$AD_pathway <- ifelse(node$node%in%ad_gene,'Yes','No')
node$AD_pathway <- factor(node$AD_pathway, levels = c('Yes','No'))
node$num_reported <- ad_gene_num$num[match(node$node, ad_gene_num$gene_name)]

ppi_select <- ppi_select %>%
  filter(gene_A %in% node$node & gene_B %in% node$node)

write.xlsx(node, 'processed_data/tmp/ppi_wmds_node.xlsx')
write.xlsx(ppi_select, 'processed_data/tmp/ppi_wmds_edge.xlsx')

# 构建图对象
ppi_graph <- tbl_graph(
  nodes = node,
  edges = ppi_select %>% rename(from = gene_A, to = gene_B),
  node_key = "node",
  directed = FALSE
)

graph <- as.igraph(ppi_graph)  # 转为 igraph 对象
layout_fr <- layout_with_fr(graph, niter = 500)
layout_fr[,1] <- layout_fr[,1] * 1.2  # x 方向缩放
layout_fr[,2] <- layout_fr[,2] * 1.2  # y 方向缩放

# 绘图
set.seed(123)
p <- ggraph(ppi_graph, layout = layout_fr) +
  geom_edge_link(color = "grey70", alpha = 0.6) +
  geom_node_point(aes(fill = AD_pathway), size = 5, shape = 21, color = "grey30") +
  geom_node_text(
    aes(label = node),
    size = 4,
    repel = TRUE,
    family = "Arial",
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.3, "lines")
  )  +
  # scale_size_continuous(range = c(2, 8)) +
  scale_fill_manual(values = c('No' = "#E0F3F8", 'Yes' = "#4DA8DA")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16, family = "Arial", face = "bold"),
    legend.text = element_text(size = 14, family = "Arial")
  )

print(p)
ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/ppi/ppi_wmds_new.png',
       p,
       width = 12,
       height = 15,
       dpi = 300)




