library(data.table)
library(readxl)
library(openxlsx)
library(writexl)
library(limma)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggprism)
library(RColorBrewer)
library(VennDiagram)
library(biomaRt)
library(patchwork)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(gground)
library(ggnewscale)
library(ggtext)
library(writexl)


source('/share/home/sunly/Rscript/cbmap_pro/cbmap_586/2_sample_process.R')
out_dir <- '/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/'


################### limma adnc
adnc_id <- intersect(sample_adnc$id, colnames(pro_50_raw)[4:ncol(pro_50_raw)])
pro_adnc_limma <- pro_50_raw[,colnames(pro_50_raw)%in%adnc_id]
pro_adnc_limma <- pro_adnc_limma[, match(adnc_id, colnames(pro_adnc_limma))]
sample_adnc <- sample_adnc[sample_adnc$id%in%adnc_id,]
sample_adnc <- sample_adnc[match(adnc_id, sample_adnc$id), ]
summary(sample_adnc$adnc_fac)

design <- model.matrix( ~ adnc_num + age + sex_male + PMD + RIN, data = sample_adnc)
fit <- lmFit(pro_adnc_limma, design)
fit <- eBayes(fit)
adnc_num_lim_re <- topTable(fit, coef = "adnc_num", number = Inf, adjust.method = "BH")
adnc_num_lim_re$gene_protein <- rownames(adnc_num_lim_re)
adnc_num_lim_re$genename <- pro_50_raw$gene_name[match(adnc_num_lim_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05
adnc_num_lim_re$protein <- pro_50_raw$protein[match(adnc_num_lim_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05

nrow(adnc_num_lim_re[adnc_num_lim_re$adj.P.Val<0.05,])
adnc_num_lim_re[adnc_num_lim_re$genename%in%c('SMOC1','CTHRC1','ICAM1','SNRPA','SNRNP70'),]
write.xlsx(adnc_num_lim_re,paste0(out_dir, 'table/dea_result/cbmap_adnc_lim_re.xlsx'))


# new report
ad_pro_meta_list <- read.xlsx('AD_pro_meta_brain_list.xlsx')
reported <- intersect(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], ad_pro_meta_list$gene[ad_pro_meta_list$FDR<0.05])
new_reported <- setdiff(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], reported)


################### limma adnc adjust cov
sample_adnc_cov <- sample_adnc
sample_adnc_cov$lbd[is.na(sample_adnc_cov$lbd)] <- 0
sample_adnc_cov$cvd[is.na(sample_adnc_cov$cvd)] <- 0
sample_adnc_cov$late[is.na(sample_adnc_cov$late)] <- 0
adnc_id_cov <- intersect(sample_adnc_cov$id, colnames(pro_50_raw)[4:ncol(pro_50_raw)])
pro_adnc_limma_cov <- pro_50_raw[,colnames(pro_50_raw)%in%adnc_id_cov]
pro_adnc_limma_cov <- pro_adnc_limma_cov[, match(adnc_id_cov, colnames(pro_adnc_limma_cov))]
sample_adnc_cov <- sample_adnc_cov[sample_adnc_cov$id%in%adnc_id_cov,]
sample_adnc_cov <- sample_adnc_cov[match(adnc_id_cov, sample_adnc_cov$id), ]
summary(sample_adnc_cov$adnc_fac)

design <- model.matrix( ~ adnc_num + age + sex_male + PMD + RIN + lbd + cvd + late, data = sample_adnc_cov)
fit <- lmFit(pro_adnc_limma_cov, design)
fit <- eBayes(fit)
adnc_num_lim_re_cov <- topTable(fit, coef = "adnc_num", number = Inf, adjust.method = "BH")
adnc_num_lim_re_cov$gene_protein <- rownames(adnc_num_lim_re_cov)
adnc_num_lim_re_cov$genename <- pro_50_raw$gene_name[match(adnc_num_lim_re_cov$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05
adnc_num_lim_re_cov$protein <- pro_50_raw$protein[match(adnc_num_lim_re_cov$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05

nrow(adnc_num_lim_re_cov[adnc_num_lim_re_cov$adj.P.Val<0.05,])
length(intersect(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], adnc_num_lim_re_cov$genename[adnc_num_lim_re_cov$adj.P.Val<0.05]))
length(setdiff(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], intersect(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], adnc_num_lim_re_cov$genename[adnc_num_lim_re_cov$adj.P.Val<0.05])))

adnc_lim_re_cov <- adnc_num_lim_re_cov
colnames(adnc_lim_re_cov)[1:5] <- paste0(colnames(adnc_lim_re_cov)[1:5],'_lbd_cvd_late')

adnc_re_cov <- adnc_num_lim_re[,c(7:9,1,5)]
adnc_re_cov<- merge(adnc_re_cov, adnc_lim_re_cov[,c(7,1,5)], by = 'gene_protein')
adnc_re_cov <- adnc_re_cov[order(adnc_re_cov$adj.P.Val),]

cor.test(adnc_re_cov$logFC, adnc_re_cov$logFC_lbd_cvd_late)

write.xlsx(adnc_re_cov, paste0(out_dir, 'table/dea_result/cbmap_adnc_lbd_cvd_lim_norm_rint.xlsx'))
write.xlsx(adnc_num_lim_re_cov, paste0(out_dir,'table/dea_result/cbmap_adnc_lim_adjlcl.xlsx'))


################### limma_braak
# braak stage(498,0-81,1-116,2-65,3-152,4-43,5-22,6-19)
braak_id <- intersect(sample_braak$id, colnames(pro_50_raw)[4:ncol(pro_50_raw)])
pro_braak_limma <- pro_50_raw[,colnames(pro_50_raw)%in%braak_id]
pro_braak_limma <- pro_braak_limma[, match(braak_id, colnames(pro_braak_limma))]
sample_braak <- sample_braak[sample_braak$id%in%braak_id,]
sample_braak <- sample_braak[match(braak_id, sample_braak$id), ]

design <- model.matrix( ~ Braak.NFT.stage + age + sex_male + PMD + RIN, data = sample_braak)
fit <- lmFit(pro_braak_limma, design)
fit <- eBayes(fit)
braak_num_lim_re <- topTable(fit, coef = "Braak.NFT.stage", number = Inf, adjust.method = "BH")
braak_num_lim_re$gene_protein <- rownames(braak_num_lim_re)
braak_num_lim_re$genename <- pro_50_raw$gene_name[match(braak_num_lim_re$gene_protein, pro_50_raw$gene_protein)] 
braak_num_lim_re$protein <- pro_50_raw$protein[match(braak_num_lim_re$gene_protein, pro_50_raw$gene_protein)] 
nrow(braak_num_lim_re[braak_num_lim_re$adj.P.Val<0.05&braak_num_lim_re$logFC<0,])

write.xlsx(braak_num_lim_re, paste0(out_dir, 'table/dea_result/cbmap_braak_lim.xlsx'))


###################### limma abeta
abeta_id <- intersect(sample_abeta$id, colnames(pro_50_raw)[4:ncol(pro_50_raw)])
pro_abeta_limma <- pro_50_raw[,colnames(pro_50_raw)%in%abeta_id]
pro_abeta_limma <- pro_abeta_limma[, match(abeta_id, colnames(pro_abeta_limma))]
sample_abeta <- sample_abeta[sample_abeta$id%in%abeta_id,]
sample_abeta <- sample_abeta[match(abeta_id, sample_abeta$id), ]

design <- model.matrix( ~ A_beta_0_3 + age + sex_male + PMD + RIN, data = sample_abeta)
fit <- lmFit(pro_abeta_limma, design)
fit <- eBayes(fit)
abeta_num_lim_re <- topTable(fit, coef = "A_beta_0_3", number = Inf, adjust.method = "BH")
abeta_num_lim_re$gene_protein <- rownames(abeta_num_lim_re)
abeta_num_lim_re$genename <- pro_50_raw$gene_name[match(abeta_num_lim_re$gene_protein, pro_50_raw$gene_protein)] 
abeta_num_lim_re$protein <- pro_50_raw$protein[match(abeta_num_lim_re$gene_protein, pro_50_raw$gene_protein)] 
nrow(abeta_num_lim_re[abeta_num_lim_re$adj.P.Val<0.05&abeta_num_lim_re$logFC>0,])

length(intersect(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05],braak_num_lim_re$genename[braak_num_lim_re$adj.P.Val<0.05]))
length(intersect(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05],abeta_num_lim_re$genename[abeta_num_lim_re$adj.P.Val<0.05]))
length(braak_num_lim_re$genename[braak_num_lim_re$adj.P.Val<0.05])
length(abeta_num_lim_re$genename[abeta_num_lim_re$adj.P.Val<0.05])
setdiff(braak_num_lim_re$genename[braak_num_lim_re$adj.P.Val<0.05],intersect(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05],braak_num_lim_re$genename[braak_num_lim_re$adj.P.Val<0.05]))
setdiff(abeta_num_lim_re$genename[abeta_num_lim_re$adj.P.Val<0.05],intersect(adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05],abeta_num_lim_re$genename[abeta_num_lim_re$adj.P.Val<0.05]))

write.xlsx(abeta_num_lim_re, paste0(out_dir, 'table/dea_result/cbmap_abeta_lim.xlsx'))


######################### enrichment 
# adnc
mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "asia")  # 或 mirror = "useast" / "asia"
gene_symbols <- unique(adnc_num_lim_re$genename)  
mapping <- getBM(attributes = c("hgnc_symbol", "entrezgene_id", "gene_biotype"),
                 filters = "hgnc_symbol",
                 values = gene_symbols,
                 mart = mart)
mapping <- mapping[!duplicated(mapping$hgnc_symbol),]

adnc_num_lim_re$entre_id <- mapping[match(adnc_num_lim_re$genename, mapping$hgnc_symbol), "entrezgene_id"]
adnc_num_lim <- adnc_num_lim_re[!is.na(adnc_num_lim_re$entre_id),]
adnc_num_lim <- adnc_num_lim[!duplicated(adnc_num_lim$entre_id),]
rank_gene_adnc_num_lim <- setNames(adnc_num_lim$logFC, adnc_num_lim$entre_id)
rank_gene_adnc_num_lim <- sort(rank_gene_adnc_num_lim, decreasing = TRUE)

categories <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
gene_sets <- do.call(rbind, lapply(categories, function(cat) {
  msigdbr(category = cat)
}))
map <- gene_sets[, c("gs_name", "entrez_gene")]
map$entrez_gene <- as.character(map$entrez_gene)

# GSAE
result_GSEA_lim = GSEA(geneList = rank_gene_adnc_num_lim, TERM2GENE = map, eps = 0)
adnc_num_gsea_lim_table <- result_GSEA_lim@result
adnc_num_gsea_lim_table$Gene_Symbols <- sapply(strsplit(adnc_num_gsea_lim_table$core_enrichment, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

adnc_num_gsea_lim_table <- as_tibble(adnc_num_gsea_lim_table)

adnc_num_gsea_lim_table <- adnc_num_gsea_lim_table %>%
  mutate(
    Direction = ifelse(NES >= 0, "Up", "Down"),
    # label_color = ifelse(NES >= 0, "red", "blue"),
    Description = gsub("_", " ", Description),
    LeadingEdgeSize = sapply(strsplit(as.character(core_enrichment), "/"), length),
    Count = LeadingEdgeSize,
    GeneRatio = LeadingEdgeSize / setSize,
    score = -log10(p.adjust) * GeneRatio,
    Description = str_wrap(Description, width = 70)
  )

format_first_upper_rest_lower <- function(x) {
  sapply(x, function(name) {
    parts <- strsplit(name, " ", fixed = TRUE)[[1]]
    if (length(parts) > 1) {
      first <- toupper(parts[1])
      rest <- tolower(paste(parts[-1], collapse = " "))
      paste(first, rest)
    } else {
      toupper(name)
    }
  })
}

adnc_num_gsea_lim_table$Description <- format_first_upper_rest_lower(adnc_num_gsea_lim_table$Description)

use_pathway <- group_by(adnc_num_gsea_lim_table, Direction) %>%
  top_n(10, wt = -p.adjust) %>%
  group_by(p.adjust) %>%
  top_n(1, wt = Count) %>%
  ungroup() %>%
  mutate(Direction = factor(Direction, levels = rev(c('Up','Down')))) %>%
  dplyr::arrange(Direction, -p.adjust) %>%
  group_by(Direction) %>%
  mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
  ungroup() %>%
  tibble::rowid_to_column('index')

xaxis_max <- max(-log10(use_pathway$pvalue)) + 1

rect.data <- group_by(use_pathway, Direction) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

pgsea <- use_pathway %>%
  ggplot(aes(-log10(pvalue), y = index, fill = Direction)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 7
  ) +
  geom_point(
    aes(x = -0.02 * xaxis_max, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.02 * xaxis_max, label = Count)
  ) +
  scale_size_continuous(name = "Count", range = c(5, 12)) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    data = data.frame(xaxis_max = xaxis_max),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  # labs(y = NULL) +
  scale_fill_manual(name = "Category", values = c("Up" = "#c2b5e3", "Down" = "#96c4a3")) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 16),
  ) +
  guides(
    size = guide_legend(override.aes = list(shape = 21)),
    fill = guide_legend(override.aes = list(shape = NA))
  )
pgsea

ggsave(paste0(out_dir, 'plot/enrichment/adnc_gsea.pdf'),
       plot = pgsea+labs(y = NULL, title = 'ADNC GSEA'),
       width = 13,
       height = 9.75)

# KEGG 
adnc_num_lim_sigpro <- adnc_num_lim$entre_id[which(adnc_num_lim$adj.P.Val < 0.05)]
kegg_result_lim <- enrichKEGG(
  gene = adnc_num_lim_sigpro,
  organism = "hsa",
  pvalueCutoff = 0.05,
  use_internal_data = FALSE
)
adnc_num_kegg_lim_table <- kegg_result_lim@result
adnc_num_kegg_lim_table$Gene_Symbols <- sapply(strsplit(adnc_num_kegg_lim_table$geneID, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

# GO
enrichgo_lim <- enrichGO(gene = adnc_num_lim_sigpro, 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "ENTREZID",  # 基因ID类型
                         ont = "ALL",
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05)
adnc_num_enrichgo_lim_table <- enrichgo_lim@result
adnc_num_enrichgo_lim_table$Gene_Symbols <- sapply(strsplit(adnc_num_enrichgo_lim_table$geneID, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

# plot
pal <- c('#c3e1e6','#f3dfb7','#dcc6dc','#96c38e')

pathway_highlight <- function(go_enrich, kegg_enrich){
  ego_readable <- setReadable(go_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  ekegg_readable <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  GO <- as.data.frame(ego_readable)
  KEGG <- as.data.frame(ekegg_readable)
  
  # 分别选取不同数量的 top 通路
  top_bp <- GO %>%
    filter(ONTOLOGY == "BP") %>%
    top_n(7, wt = -pvalue)
  
  top_cc <- GO %>%
    filter(ONTOLOGY == "CC") %>%
    top_n(5, wt = -pvalue)
  
  top_mf <- GO %>%
    filter(ONTOLOGY == "MF") %>%
    top_n(6, wt = -pvalue)
  
  top_kegg <- KEGG %>%
    top_n(2, wt = -pvalue) %>%
    mutate(ONTOLOGY = "KEGG")
  
  # 合并所有
  use_pathway <- bind_rows(top_bp, top_cc, top_mf, top_kegg) %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
    arrange(ONTOLOGY, -pvalue) %>%
    group_by(ONTOLOGY) %>%
    mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
    tibble::rowid_to_column("index")
  
  xaxis_max <- max(-log10(use_pathway$pvalue)) + 1
  
  rect.data <- group_by(use_pathway, ONTOLOGY) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(
      ymax = cumsum(n),
      ymin = lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )
  
  p <- use_pathway %>%
    ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Description), width = 0.6, alpha = 0.8
    ) +
    geom_text(
      aes(x = 0.05, label = Description),
      hjust = 0, size = 7
    ) +
    geom_point(
      aes(x = -0.02 * xaxis_max, size = Count),
      shape = 21
    ) +
    geom_text(
      aes(x = -0.02 * xaxis_max, label = Count)
    ) +
    scale_size_continuous(name = "Count", range = c(5, 12)) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      data = data.frame(xaxis_max = xaxis_max),
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(name = "Category", values = pal) +
    scale_colour_manual(values = pal) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(),
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 16)
    )+
    guides(
      size = guide_legend(override.aes = list(shape = 21)),
      fill = guide_legend(override.aes = list(shape = NA))
    )
  print(p)
  
  return(p=p)
}

pgk <- pathway_highlight(enrichgo_lim, kegg_result_lim)

ggsave(paste0(out_dir, 'plot/enrichment/adnc_go_kegg_ga.pdf'),
       plot = pgk+labs(y = NULL, title = 'ADNC GO & KEGG'),
       width = 12,
       height = 9)

enrichment_result <- list(
  adnc_num_gsea_lim = adnc_num_gsea_lim_table,
  adnc_num_kegg_lim = adnc_num_kegg_lim_table,
  adnc_num_enrichgo_lim = adnc_num_enrichgo_lim_table
)

write_xlsx(enrichment_result, paste0(out_dir, 'table/enrich_result/cbmap_adnc_gkg.xlsx'))

# abeta
mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "asia")  # 或 mirror = "useast" / "asia"
gene_symbols <- unique(abeta_num_lim_re$genename)  
mapping <- getBM(attributes = c("hgnc_symbol", "entrezgene_id", "gene_biotype"),
                 filters = "hgnc_symbol",
                 values = gene_symbols,
                 mart = mart)
mapping <- mapping[!duplicated(mapping$hgnc_symbol),]

abeta_num_lim_re$entre_id <- mapping[match(abeta_num_lim_re$genename, mapping$hgnc_symbol), "entrezgene_id"]
abeta_num_lim <- abeta_num_lim_re[!is.na(abeta_num_lim_re$entre_id),]
abeta_num_lim <- abeta_num_lim[!duplicated(abeta_num_lim$entre_id),]
rank_gene_abeta_num_lim <- setNames(abeta_num_lim$logFC, abeta_num_lim$entre_id)
rank_gene_abeta_num_lim <- sort(rank_gene_abeta_num_lim, decreasing = TRUE)

categories <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
gene_sets <- do.call(rbind, lapply(categories, function(cat) {
  msigdbr(category = cat)
}))
map <- gene_sets[, c("gs_name", "entrez_gene")]
map$entrez_gene <- as.character(map$entrez_gene)

# GSAE
result_GSEA_lim = GSEA(geneList = rank_gene_abeta_num_lim, TERM2GENE = map, eps = 0)
abeta_num_gsea_lim_table <- result_GSEA_lim@result
abeta_num_gsea_lim_table$Gene_Symbols <- sapply(strsplit(abeta_num_gsea_lim_table$core_enrichment, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

abeta_num_gsea_lim_table <- as_tibble(abeta_num_gsea_lim_table)

abeta_num_gsea_lim_table <- abeta_num_gsea_lim_table %>%
  mutate(
    Direction = ifelse(NES >= 0, "Up", "Down"),
    # label_color = ifelse(NES >= 0, "red", "blue"),
    Description = gsub("_", " ", Description),
    LeadingEdgeSize = sapply(strsplit(as.character(core_enrichment), "/"), length),
    Count = LeadingEdgeSize,
    GeneRatio = LeadingEdgeSize / setSize,
    score = -log10(p.adjust) * GeneRatio,
    Description = str_wrap(Description, width = 70)
  )

format_first_upper_rest_lower <- function(x) {
  sapply(x, function(name) {
    parts <- strsplit(name, " ", fixed = TRUE)[[1]]
    if (length(parts) > 1) {
      first <- toupper(parts[1])
      rest <- tolower(paste(parts[-1], collapse = " "))
      paste(first, rest)
    } else {
      toupper(name)
    }
  })
}

abeta_num_gsea_lim_table$Description <- format_first_upper_rest_lower(abeta_num_gsea_lim_table$Description)

use_pathway <- group_by(abeta_num_gsea_lim_table, Direction) %>%
  top_n(10, wt = -p.adjust) %>%
  group_by(p.adjust) %>%
  top_n(1, wt = Count) %>%
  ungroup() %>%
  mutate(Direction = factor(Direction, levels = rev(c('Up','Down')))) %>%
  dplyr::arrange(Direction, -p.adjust) %>%
  group_by(Direction) %>%
  mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
  ungroup() %>%
  tibble::rowid_to_column('index')

xaxis_max <- max(-log10(use_pathway$pvalue)) + 1

rect.data <- group_by(use_pathway, Direction) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

pgsea <- use_pathway %>%
  ggplot(aes(-log10(pvalue), y = index, fill = Direction)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 7
  ) +
  geom_point(
    aes(x = -0.02 * xaxis_max, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.02 * xaxis_max, label = Count)
  ) +
  scale_size_continuous(name = "Count", range = c(5, 12)) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    data = data.frame(xaxis_max = xaxis_max),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  # labs(y = NULL) +
  scale_fill_manual(name = "Category", values = c("Up" = "#c2b5e3", "Down" = "#96c4a3")) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 16),
  ) +
  guides(
    size = guide_legend(override.aes = list(shape = 21)),
    fill = guide_legend(override.aes = list(shape = NA))
  )
pgsea

ggsave(paste0(out_dir, 'plot/enrichment/abeta_gsea.pdf'),
       plot = pgsea+labs(y = NULL, title = 'Abeta GSEA'),
       width = 13,
       height = 9.75)

# KEGG 
abeta_num_lim_sigpro <- abeta_num_lim$entre_id[which(abeta_num_lim$adj.P.Val < 0.05)]
kegg_result_lim <- enrichKEGG(
  gene = abeta_num_lim_sigpro,
  organism = "hsa",
  pvalueCutoff = 0.05,
  use_internal_data = FALSE
)
abeta_num_kegg_lim_table <- kegg_result_lim@result
abeta_num_kegg_lim_table$Gene_Symbols <- sapply(strsplit(abeta_num_kegg_lim_table$geneID, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

# GO
enrichgo_lim <- enrichGO(gene = abeta_num_lim_sigpro, 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "ENTREZID",  # 基因ID类型
                         ont = "ALL",
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05)
abeta_num_enrichgo_lim_table <- enrichgo_lim@result
abeta_num_enrichgo_lim_table$Gene_Symbols <- sapply(strsplit(abeta_num_enrichgo_lim_table$geneID, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

# plot
pal <- c('#c3e1e6','#f3dfb7','#dcc6dc','#96c38e')

pathway_highlight <- function(go_enrich, kegg_enrich){
  ego_readable <- setReadable(go_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  ekegg_readable <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  GO <- as.data.frame(ego_readable)
  KEGG <- as.data.frame(ekegg_readable)
  
  # 分别选取不同数量的 top 通路
  top_bp <- GO %>%
    filter(ONTOLOGY == "BP") %>%
    top_n(6, wt = -pvalue)
  
  top_cc <- GO %>%
    filter(ONTOLOGY == "CC") %>%
    top_n(5, wt = -pvalue)
  
  top_mf <- GO %>%
    filter(ONTOLOGY == "MF") %>%
    top_n(5, wt = -pvalue)
  
  top_kegg <- KEGG %>%
    top_n(4, wt = -pvalue) %>%
    mutate(ONTOLOGY = "KEGG")
  
  # 合并所有
  use_pathway <- bind_rows(top_bp, top_cc, top_mf, top_kegg) %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
    arrange(ONTOLOGY, -pvalue) %>%
    group_by(ONTOLOGY) %>%
    mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
    tibble::rowid_to_column("index")
  
  xaxis_max <- max(-log10(use_pathway$pvalue)) + 1
  
  rect.data <- group_by(use_pathway, ONTOLOGY) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(
      ymax = cumsum(n),
      ymin = lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )
  
  p <- use_pathway %>%
    ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Description), width = 0.6, alpha = 0.8
    ) +
    geom_text(
      aes(x = 0.05, label = Description),
      hjust = 0, size = 7
    ) +
    geom_point(
      aes(x = -0.02 * xaxis_max, size = Count),
      shape = 21
    ) +
    geom_text(
      aes(x = -0.02 * xaxis_max, label = Count)
    ) +
    scale_size_continuous(name = "Count", range = c(5, 12)) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      data = data.frame(xaxis_max = xaxis_max),
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(name = "Category", values = pal) +
    scale_colour_manual(values = pal) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(),
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 16)
    )+
    guides(
      size = guide_legend(override.aes = list(shape = 21)),
      fill = guide_legend(override.aes = list(shape = NA))
    )
  print(p)
  
  return(p=p)
}

pgk <- pathway_highlight(enrichgo_lim, kegg_result_lim)

ggsave(paste0(out_dir, 'plot/enrichment/abeta_go_kegg_ga.pdf'),
       plot = pgk+labs(y = NULL, title = 'Abeta GO & KEGG'),
       width = 12,
       height = 9)

enrichment_result <- list(
  abeta_num_gsea_lim = abeta_num_gsea_lim_table,
  abeta_num_kegg_lim = abeta_num_kegg_lim_table,
  abeta_num_enrichgo_lim = abeta_num_enrichgo_lim_table
)

write_xlsx(enrichment_result, paste0(out_dir, 'table/enrich_result/cbmap_abeta_gkg.xlsx'))

# braak
mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "asia")  # 或 mirror = "useast" / "asia"
gene_symbols <- unique(braak_num_lim_re$genename)  
mapping <- getBM(attributes = c("hgnc_symbol", "entrezgene_id", "gene_biotype"),
                 filters = "hgnc_symbol",
                 values = gene_symbols,
                 mart = mart)
mapping <- mapping[!duplicated(mapping$hgnc_symbol),]

braak_num_lim_re$entre_id <- mapping[match(braak_num_lim_re$genename, mapping$hgnc_symbol), "entrezgene_id"]
braak_num_lim <- braak_num_lim_re[!is.na(braak_num_lim_re$entre_id),]
braak_num_lim <- braak_num_lim[!duplicated(braak_num_lim$entre_id),]
rank_gene_braak_num_lim <- setNames(braak_num_lim$logFC, braak_num_lim$entre_id)
rank_gene_braak_num_lim <- sort(rank_gene_braak_num_lim, decreasing = TRUE)

categories <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
gene_sets <- do.call(rbind, lapply(categories, function(cat) {
  msigdbr(category = cat)
}))
map <- gene_sets[, c("gs_name", "entrez_gene")]
map$entrez_gene <- as.character(map$entrez_gene)

# GSAE
result_GSEA_lim = GSEA(geneList = rank_gene_braak_num_lim, TERM2GENE = map, eps = 0)
braak_num_gsea_lim_table <- result_GSEA_lim@result
braak_num_gsea_lim_table$Gene_Symbols <- sapply(strsplit(braak_num_gsea_lim_table$core_enrichment, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

braak_num_gsea_lim_table <- as_tibble(braak_num_gsea_lim_table)

braak_num_gsea_lim_table <- braak_num_gsea_lim_table %>%
  mutate(
    Direction = ifelse(NES >= 0, "Up", "Down"),
    # label_color = ifelse(NES >= 0, "red", "blue"),
    Description = gsub("_", " ", Description),
    LeadingEdgeSize = sapply(strsplit(as.character(core_enrichment), "/"), length),
    Count = LeadingEdgeSize,
    GeneRatio = LeadingEdgeSize / setSize,
    score = -log10(p.adjust) * GeneRatio,
    Description = str_wrap(Description, width = 70)
  )

format_first_upper_rest_lower <- function(x) {
  sapply(x, function(name) {
    parts <- strsplit(name, " ", fixed = TRUE)[[1]]
    if (length(parts) > 1) {
      first <- toupper(parts[1])
      rest <- tolower(paste(parts[-1], collapse = " "))
      paste(first, rest)
    } else {
      toupper(name)
    }
  })
}

braak_num_gsea_lim_table$Description <- format_first_upper_rest_lower(braak_num_gsea_lim_table$Description)

use_pathway <- group_by(braak_num_gsea_lim_table, Direction) %>%
  top_n(10, wt = -p.adjust) %>%
  group_by(p.adjust) %>%
  top_n(1, wt = Count) %>%
  ungroup() %>%
  mutate(Direction = factor(Direction, levels = rev(c('Up','Down')))) %>%
  dplyr::arrange(Direction, -p.adjust) %>%
  group_by(Direction) %>%
  mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
  ungroup() %>%
  tibble::rowid_to_column('index')

xaxis_max <- max(-log10(use_pathway$pvalue)) + 1

rect.data <- group_by(use_pathway, Direction) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

pgsea <- use_pathway %>%
  ggplot(aes(-log10(pvalue), y = index, fill = Direction)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 7
  ) +
  geom_point(
    aes(x = -0.02 * xaxis_max, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.02 * xaxis_max, label = Count)
  ) +
  scale_size_continuous(name = "Count", range = c(5, 12)) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    data = data.frame(xaxis_max = xaxis_max),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  # labs(y = NULL) +
  scale_fill_manual(name = "Category", values = c("Up" = "#c2b5e3", "Down" = "#96c4a3")) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 16),
  ) +
  guides(
    size = guide_legend(override.aes = list(shape = 21)),
    fill = guide_legend(override.aes = list(shape = NA))
  )
pgsea

ggsave(paste0(out_dir, 'plot/enrichment/braak_gsea.pdf'),
       plot = pgsea+labs(y = NULL, title = 'braak GSEA'),
       width = 13,
       height = 9.75)

# KEGG 
braak_num_lim_sigpro <- braak_num_lim$entre_id[which(braak_num_lim$adj.P.Val < 0.05)]
kegg_result_lim <- enrichKEGG(
  gene = braak_num_lim_sigpro,
  organism = "hsa",
  pvalueCutoff = 0.05,
  use_internal_data = FALSE
)
braak_num_kegg_lim_table <- kegg_result_lim@result
braak_num_kegg_lim_table$Gene_Symbols <- sapply(strsplit(braak_num_kegg_lim_table$geneID, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

# GO
enrichgo_lim <- enrichGO(gene = braak_num_lim_sigpro, 
                         OrgDb = org.Hs.eg.db, 
                         keyType = "ENTREZID",  # 基因ID类型
                         ont = "ALL",
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 0.05)
braak_num_enrichgo_lim_table <- enrichgo_lim@result
braak_num_enrichgo_lim_table$Gene_Symbols <- sapply(strsplit(braak_num_enrichgo_lim_table$geneID, "/"), function(x) {
  symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
  paste(symbols, collapse = ", ")
})

# plot
pal <- c('#c3e1e6','#f3dfb7','#dcc6dc','#96c38e')

pathway_highlight <- function(go_enrich, kegg_enrich){
  ego_readable <- setReadable(go_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  ekegg_readable <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  GO <- as.data.frame(ego_readable)
  KEGG <- as.data.frame(ekegg_readable)
  
  # 分别选取不同数量的 top 通路
  top_bp <- GO %>%
    filter(ONTOLOGY == "BP") %>%
    top_n(5, wt = -pvalue)
  
  top_cc <- GO %>%
    filter(ONTOLOGY == "CC") %>%
    top_n(5, wt = -pvalue)
  
  top_mf <- GO %>%
    filter(ONTOLOGY == "MF") %>%
    top_n(5, wt = -pvalue)
  
  top_kegg <- KEGG %>%
    top_n(5, wt = -pvalue) %>%
    mutate(ONTOLOGY = "KEGG")
  
  # 合并所有
  use_pathway <- bind_rows(top_bp, top_cc, top_mf, top_kegg) %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
    arrange(ONTOLOGY, -pvalue) %>%
    group_by(ONTOLOGY) %>%
    mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
    tibble::rowid_to_column("index")
  
  xaxis_max <- max(-log10(use_pathway$pvalue)) + 1
  
  rect.data <- group_by(use_pathway, ONTOLOGY) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(
      ymax = cumsum(n),
      ymin = lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )
  
  p <- use_pathway %>%
    ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
    geom_round_col(
      aes(y = Description), width = 0.6, alpha = 0.8
    ) +
    geom_text(
      aes(x = 0.05, label = Description),
      hjust = 0, size = 7
    ) +
    geom_point(
      aes(x = -0.02 * xaxis_max, size = Count),
      shape = 21
    ) +
    geom_text(
      aes(x = -0.02 * xaxis_max, label = Count)
    ) +
    scale_size_continuous(name = "Count", range = c(5, 12)) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      data = data.frame(xaxis_max = xaxis_max),
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(name = "Category", values = pal) +
    scale_colour_manual(values = pal) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(),
      axis.text.x = element_text(size = 15),
      axis.title.x = element_text(size = 16)
    )+
    guides(
      size = guide_legend(override.aes = list(shape = 21)),
      fill = guide_legend(override.aes = list(shape = NA))
    )
  print(p)
  
  return(p=p)
}

pgk <- pathway_highlight(enrichgo_lim, kegg_result_lim)

ggsave(paste0(out_dir, 'plot/enrichment/braak_go_kegg_ga.pdf'),
       plot = pgk+labs(y = NULL, title = 'braak GO & KEGG'),
       width = 12,
       height = 9)

enrichment_result <- list(
  braak_num_gsea_lim = braak_num_gsea_lim_table,
  braak_num_kegg_lim = braak_num_kegg_lim_table,
  braak_num_enrichgo_lim = braak_num_enrichgo_lim_table
)

write_xlsx(enrichment_result, paste0(out_dir, 'table/enrich_result/cbmap_braak_gkg.xlsx'))


########################## sensitivity analysis
# adnc
sample_adnc_sen <- sample_adnc[sample_adnc$bank%in%c('zju','pumc'),]
adnc_id_sen <- intersect(colnames(pro_50_raw),sample_adnc_sen$id)

pro_adnc_sen_limma <- pro_50_raw[,colnames(pro_50_raw)%in%adnc_id_sen]
pro_adnc_sen_limma <- pro_adnc_sen_limma[, match(adnc_id_sen, colnames(pro_adnc_sen_limma))]
sample_adnc_sen <- sample_adnc_sen[sample_adnc_sen$id%in%adnc_id_sen,]
sample_adnc_sen <- sample_adnc_sen[match(adnc_id_sen, sample_adnc_sen$id), ]
summary(sample_adnc_sen$adnc_fac)

design <- model.matrix( ~ adnc_num + age + sex_male + PMD + RIN, data = sample_adnc_sen)
fit <- lmFit(pro_adnc_sen_limma, design)
fit <- eBayes(fit)
adnc_num_lim_sen_re <- topTable(fit, coef = "adnc_num", number = Inf, adjust.method = "BH")
adnc_num_lim_sen_re$gene_protein <- rownames(adnc_num_lim_sen_re)
adnc_num_lim_sen_re$genename <- pro_50_raw$gene_name[match(adnc_num_lim_sen_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05
adnc_num_lim_sen_re$protein <- pro_50_raw$protein[match(adnc_num_lim_sen_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05

nrow(adnc_num_lim_sen_re[adnc_num_lim_sen_re$adj.P.Val<0.05,])
adnc_num_lim_sen_re[adnc_num_lim_sen_re$genename%in%c('SMOC1','CTHRC1','ICAM1','SNRPA','SNRNP70'),]
write.xlsx(adnc_num_lim_sen_re,paste0(out_dir, 'table/dea_result/cbmap_adnc_lim_sen_re.xlsx'))

# braak
sample_braak_sen <- sample_braak[sample_braak$bank%in%c('zju','pumc'),]
braak_id_sen <- intersect(colnames(pro_50_raw),sample_braak_sen$id)

pro_braak_sen_limma <- pro_50_raw[,colnames(pro_50_raw)%in%braak_id_sen]
pro_braak_sen_limma <- pro_braak_sen_limma[, match(braak_id_sen, colnames(pro_braak_sen_limma))]
sample_braak_sen <- sample_braak_sen[sample_braak_sen$id%in%braak_id_sen,]
sample_braak_sen <- sample_braak_sen[match(braak_id_sen, sample_braak_sen$id), ]
summary(sample_braak_sen$braak_fac)

design <- model.matrix( ~ Braak.NFT.stage + age + sex_male + PMD + RIN, data = sample_braak_sen)
fit <- lmFit(pro_braak_sen_limma, design)
fit <- eBayes(fit)
braak_num_lim_sen_re <- topTable(fit, coef = "Braak.NFT.stage", number = Inf, adjust.method = "BH")
braak_num_lim_sen_re$gene_protein <- rownames(braak_num_lim_sen_re)
braak_num_lim_sen_re$genename <- pro_50_raw$gene_name[match(braak_num_lim_sen_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05
braak_num_lim_sen_re$protein <- pro_50_raw$protein[match(braak_num_lim_sen_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05

nrow(braak_num_lim_sen_re[braak_num_lim_sen_re$adj.P.Val<0.05,])
braak_num_lim_sen_re[braak_num_lim_sen_re$genename%in%c('SMOC1','CTHRC1','ICAM1','SNRPA','SNRNP70'),]
write.xlsx(braak_num_lim_sen_re,paste0(out_dir, 'table/dea_result/cbmap_braak_lim_sen_re.xlsx'))

# abeta
sample_abeta_sen <- sample_abeta[sample_abeta$bank%in%c('zju','pumc'),]
abeta_id_sen <- intersect(colnames(pro_50_raw),sample_abeta_sen$id)

pro_abeta_sen_limma <- pro_50_raw[,colnames(pro_50_raw)%in%abeta_id_sen]
pro_abeta_sen_limma <- pro_abeta_sen_limma[, match(abeta_id_sen, colnames(pro_abeta_sen_limma))]
sample_abeta_sen <- sample_abeta_sen[sample_abeta_sen$id%in%abeta_id_sen,]
sample_abeta_sen <- sample_abeta_sen[match(abeta_id_sen, sample_abeta_sen$id), ]
summary(sample_abeta_sen$abeta_fac)

design <- model.matrix( ~ A_beta_0_3 + age + sex_male + PMD + RIN, data = sample_abeta_sen)
fit <- lmFit(pro_abeta_sen_limma, design)
fit <- eBayes(fit)
abeta_num_lim_sen_re <- topTable(fit, coef = "A_beta_0_3", number = Inf, adjust.method = "BH")
abeta_num_lim_sen_re$gene_protein <- rownames(abeta_num_lim_sen_re)
abeta_num_lim_sen_re$genename <- pro_50_raw$gene_name[match(abeta_num_lim_sen_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05
abeta_num_lim_sen_re$protein <- pro_50_raw$protein[match(abeta_num_lim_sen_re$gene_protein, pro_50_raw$gene_protein)] # PLCB1过fdr校正，P值均小于0.05

nrow(abeta_num_lim_sen_re[abeta_num_lim_sen_re$adj.P.Val<0.05,])
abeta_num_lim_sen_re[abeta_num_lim_sen_re$genename%in%c('SMOC1','CTHRC1','ICAM1','SNRPA','SNRNP70'),]
write.xlsx(abeta_num_lim_sen_re,paste0(out_dir, 'table/dea_result/cbmap_abeta_lim_sen_re.xlsx'))


########################## miami plot
pathway_func <- read.xlsx('splicing_pathway_categories.xlsx')
sheet_names <- excel_sheets(paste0(out_dir,'table/enrich_result/cbmap_adnc_gkg.xlsx'))
enrichment <- lapply(sheet_names, function(sheet) {
  read_excel(paste0(out_dir,'table/enrich_result/cbmap_adnc_gkg.xlsx'),sheet = sheet)
})

adnc_go <- enrichment[[3]]
adnc_go$Category <- pathway_func$Category[match(adnc_go$Description, pathway_func$Term)]

category <- unique(adnc_go$Category)
gene_categray <- adnc_go %>%
  separate_rows(geneID, sep = "/") %>%
  distinct(Category, geneID)  # 去重
gene_categray$Category <- ifelse(gene_categray$Category=="Cytoskeleton & contractile structure", "Cytoskeleton", gene_categray$Category)
gene_categray$cat <- ifelse(gene_categray$Category%in%c("RNA processing & splicing","Synapse & neuron","Cytoskeleton","Cell adhesion & signaling"), gene_categray$Category, NA)

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = 'www')
genes <- unique(adnc_num_lim_re$genename)
gene_anno <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = genes,
  mart = ensembl
)

valid_chr <- c(as.character(1:22), "X", "Y")
gene_anno <- gene_anno[gene_anno$chromosome_name %in% valid_chr, ]
gene_anno$chromosome_name <- factor(gene_anno$chromosome_name, levels = valid_chr)

dual_miami_plot <- function(data, gene_anno, fdr_col = "P.Value", logfc_col = "logFC", genename_col = "genename",
                            chr_col = "chromosome_name", start_col = "start_position", end_col = "end_position",
                            top_n = 30) {
  
  fdr_thresh =  max(data$P.Value[data$adj.P.Val<0.05])
  
  # 合并注释
  data <- merge(data, gene_anno, by.x = genename_col, by.y = "hgnc_symbol")
  
  # 保留有效染色体
  valid_chr <- c(as.character(1:22), "X", "Y")
  data <- data[data[[chr_col]] %in% valid_chr, ]
  data[[chr_col]] <- factor(as.character(data[[chr_col]]), levels = valid_chr)
  
  # 染色体颜色
  my_upper_colors <- c('grey40','grey')
  my_lower_colors <- c('grey40','grey')
  
  chr_colors_up <- rep(my_upper_colors, length.out = length(valid_chr))
  names(chr_colors_up) <- valid_chr
  chr_colors_down <- rep(my_lower_colors, length.out = length(valid_chr))
  names(chr_colors_down) <- valid_chr
  
  # 构建完整染色体信息表，即使部分染色体无数据也保留
  all_chr_info <- gene_anno %>%
    filter(chromosome_name %in% valid_chr) %>%
    group_by(chromosome_name) %>%
    summarise(chr_len = max(end_position, na.rm = TRUE)) %>%
    ungroup() %>%
    complete(chromosome_name = factor(valid_chr, levels = valid_chr), fill = list(chr_len = 1)) %>%
    mutate(chr_start = cumsum(as.numeric(lag(chr_len, default = 0))))
  
  # breaks 与 labels
  x_breaks <- all_chr_info$chr_start + all_chr_info$chr_len / 2
  x_labels <- all_chr_info$chromosome_name
  x_limits <- c(0, sum(all_chr_info$chr_len, na.rm = TRUE))
  
  # 加入 pos_cum 与 logp
  colnames(all_chr_info)[colnames(all_chr_info) == "chromosome_name"] <- chr_col
  
  data <- left_join(data, all_chr_info, by = chr_col) %>%
    mutate(
      pos_cum = as.numeric(.data[[start_col]]) + chr_start,
      logp = -log10(.data[[fdr_col]]),
      logp_dir = ifelse(.data[[logfc_col]] >= 0, logp, -logp)
    )
  
  # 加入富集功能
  gene_cat <- gene_categray[!is.na(gene_categray$cat),]
  
  priority <- c("RNA processing & splicing", "Cytoskeleton", "Synapse & neuron", "Cell adhesion & signaling")
  gene_cat <- gene_cat %>%
    dplyr::mutate(Category = factor(Category, levels = priority)) %>%
    dplyr::arrange(geneID, Category) %>%
    dplyr::group_by(geneID) %>%
    dplyr::slice(1)
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
  res <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_cat$geneID,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
  )
  
  gene_cat$genename <- res$SYMBOL[match(gene_cat$geneID, res$ENTREZID)]
  
  data <- data %>%
    left_join(gene_cat[,c(1,4)], by = c("genename" = "genename"))
  data <- data[order(data$adj.P.Val),]
  top_pro <- data$genename[1:top_n]
  
  # Top hits 和阈值
  top_hits_pos <- data %>%
    filter(.data[[logfc_col]] >= 0) %>%
    arrange(.data[[fdr_col]]) %>%
    head(10)
  
  top_hits_neg <- data %>%
    filter(.data[[logfc_col]] < 0) %>%
    arrange(.data[[fdr_col]]) %>%
    head(10)
  
  top_hits <- bind_rows(top_hits_pos, top_hits_neg)
  
  select_data <- data[data$genename%in%top_pro & !is.na(data$Category),]
  top_hits <- rbind(top_hits, select_data)
  top_hits <- unique(top_hits)
  top_hits <- top_hits %>%
    mutate(highlight_fill = ifelse(is.na(Category), NA, as.character(Category)))
  
  threshold <- -log10(fdr_thresh)
  maxp <- max(abs(data$logp_dir), na.rm = TRUE)
  
  categories <- na.omit(unique(top_hits$Category))
  category_colors <- c(
    "RNA processing & splicing" = "#C3D9F0",        
    "Cytoskeleton" = "#FFF2C2", 
    "Synapse & neuron" = "#B5E6DE",                 
    "Cell adhesion & signaling" = "#D3C9EB"         
  )
  
  category_colors <- category_colors[categories]
  
  # upper plot
  upper_plot <- ggplot(data = data[data$logp_dir > 0, ],
                       aes(x = pos_cum, y = logp_dir, color = .data[[chr_col]])) +
    # 背景点（按染色体着色）
    geom_point(alpha = 0.7, size = 1, show.legend = FALSE) +
    scale_color_manual(values = chr_colors_up) +
    
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed") +
    
    # 标签1：有分类，颜色按 Category
    geom_label_repel(
      data = subset(top_hits, logp_dir > 0 & !is.na(Category)),
      aes(x = pos_cum, y = logp_dir, label = .data[[genename_col]], fill = Category),
      color = 'black',
      size = 5,
      segment.size = 0.2,
      point.padding = 0.7,
      box.padding = 0.7,
      max.overlaps = Inf,
      force = 10,
      force_pull = 0.05,
      min.segment.length = 0,
      segment.curvature = 0.1,
      segment.ncp = 3,
      direction = "both"
      # show.legend = TRUE
    ) +
    
    # 标签2：无分类，颜色按染色体
    geom_label_repel(
      data = subset(top_hits, logp_dir > 0 & is.na(Category)),
      aes(x = pos_cum, y = logp_dir, label = .data[[genename_col]], fill = 'white'),
      color = 'black',
      size = 5,
      segment.size = 0.2,
      point.padding = 0.7,
      box.padding = 0.7,
      max.overlaps = Inf,
      force = 10,
      force_pull = 0.05,
      min.segment.length = 0,
      segment.curvature = 0.1,
      segment.ncp = 3,
      direction = "both"
      # show.legend = TRUE
    ) +
    
    # 图例只显示 Category，不显示染色体颜色图例
    scale_fill_manual(values = category_colors, na.translate = FALSE, name = NULL, drop = FALSE, guide = 'none') +
    
    scale_y_continuous(name = "Positive -log10(P)", limits = c(0, maxp)) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    coord_cartesian(xlim = x_limits) +
    
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(b = 0, l = 10),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18)
    )
  
  # lower plot
  lower_plot <- ggplot(data = data[data$logp_dir < 0, ],
                       aes(x = pos_cum, y = logp_dir, color = .data[[chr_col]])) +
    # 背景点（按染色体着色）
    geom_point(alpha = 0.7, size = 1, show.legend = FALSE) +
    scale_color_manual(values = chr_colors_down) +
    
    geom_hline(yintercept = -threshold, color = "red", linetype = "dashed") +
    
    # 标签1：有分类，颜色按 Category
    geom_label_repel(
      data = subset(top_hits, logp_dir < 0 & !is.na(Category)),
      aes(x = pos_cum, y = logp_dir, label = .data[[genename_col]], fill = Category),
      color = 'black',
      size = 5,
      segment.size = 0.2,
      point.padding = 0.7,
      box.padding = 0.7,
      max.overlaps = Inf,
      force = 10,
      force_pull = 0.05,
      min.segment.length = 0,
      segment.curvature = 0.1,
      segment.ncp = 3,
      direction = "both"
      # show.legend = TRUE
    ) +
    
    # 标签2：无分类，颜色按染色体
    geom_label_repel(
      data = subset(top_hits, logp_dir < 0 & is.na(Category)),
      aes(x = pos_cum, y = logp_dir, label = .data[[genename_col]], fill = 'white'),
      color = 'black',
      size = 5,
      segment.size = 0.2,
      point.padding = 0.7,
      box.padding = 0.7,
      max.overlaps = Inf,
      force = 10,
      force_pull = 0.05,
      min.segment.length = 0,
      segment.curvature = 0.1,
      segment.ncp = 3,
      direction = "both"
      # show.legend = TRUE
    ) +
    
    #图例只显示 Category，不显示染色体颜色图例
    scale_fill_manual(values = category_colors,
                      na.translate = FALSE,
                      name = NULL,
                      drop = FALSE,
                      guide = guide_legend( # 同样修改
                        override.aes = list(label = "", size = 0),
                        keywidth = unit(12, "pt"),
                        keyheight = unit(12, "pt")
                      )) +
    
    scale_y_continuous(name = "Negative -log10(P)", limits = c(-maxp,0),labels = abs) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, position = 'top') +
    coord_cartesian(xlim = x_limits) +
    
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      plot.margin = margin(t = 0, l = 10),
      # plot.margin = margin(t = 5, r = 5, b = 5, l = 10),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      legend.text = element_text(size = 13)
    )
  
  # 合并图形
  combined_plot <- upper_plot / lower_plot + plot_layout(guides = "collect")
  combined_plot <- combined_plot & theme(legend.position = "bottom", 
                                         legend.box = "horizontal")
  
  print(combined_plot)
  return(combined_plot)
}

adnc_plot <- dual_miami_plot(adnc_num_lim_re, gene_anno)
ggsave(paste0(out_dir, 'plot/miami/adnc_protein.pdf'),
       adnc_plot,
       height = 9,
       width = 12)

adnc_cov_plot <- dual_miami_plot(adnc_num_lim_re_cov, gene_anno)
ggsave(paste0(out_dir, 'plot/miami/adnc_protein_cov.pdf'),
       adnc_cov_plot,
       height = 9,
       width = 12)
adnc_sen_plot <- dual_miami_plot(adnc_num_lim_sen_re, gene_anno)
ggsave(paste0(out_dir, 'plot/miami/adnc_protein_sen.pdf'),
       adnc_sen_plot,
       height = 9,
       width = 12)

braak_plot <- dual_miami_plot(braak_num_lim_re, gene_anno)
ggsave(paste0(out_dir, 'plot/miami/braak_protein.pdf'),
       braak_plot,
       height = 9,
       width = 12) 
braak_sen_plot <- dual_miami_plot(braak_num_lim_sen_re, gene_anno)
ggsave(paste0(out_dir, 'plot/miami/braak_protein_sen.pdf'),
       braak_sen_plot,
       height = 9,
       width = 12) 
abeta_plot <- dual_miami_plot(abeta_num_lim_re, gene_anno)
ggsave(paste0(out_dir, 'plot/miami/abeta_protein.pdf'),
       abeta_plot,
       height = 9,
       width = 12) 
abeta_sen_plot <- dual_miami_plot(abeta_num_lim_sen_re, gene_anno)
ggsave(paste0(out_dir, 'plot/miami/abeta_protein_sen.pdf'),
       abeta_sen_plot,
       height = 9,
       width = 12) 


######################### top gene box plot
pro_50_log <- as.data.frame(fread('processed_data/human_brain_pro_50_log_729.txt'))
rownames(pro_50_log) <- pro_50_log$protein

# adnc
pro_50_log_adnc <- pro_50_log[,colnames(pro_50_log)%in%sample_adnc$id]
pro_50_log_adnc_t <- data.frame(t(pro_50_log_adnc))
pro_50_log_adnc_t$id <- rownames(pro_50_log_adnc_t)
pro_50_log_adnc_t$adnc <- sample_adnc$adnc_num[match(pro_50_log_adnc_t$id, sample_adnc$id)]
pro_50_log_adnc_t$bank <- sample_adnc$bank[match(pro_50_log_adnc_t$id, sample_adnc$id)]
pro_50_log_adnc_t$adnc_type <- ifelse(pro_50_log_adnc_t$adnc==0,'HC',
                                      ifelse(pro_50_log_adnc_t$adnc==1,'ADNC-L',
                                             ifelse(pro_50_log_adnc_t$adnc==2,'ADNC-M','ADNC-H')))
pro_50_log_adnc_t$adnc_type <- factor(pro_50_log_adnc_t$adnc_type, levels = c('HC','ADNC-L','ADNC-M','ADNC-H'))

# braak
pro_50_log_braak <- pro_50_log[,colnames(pro_50_log)%in%sample_braak$id]
pro_50_log_braak_t <- data.frame(t(pro_50_log_braak))
pro_50_log_braak_t$id <- rownames(pro_50_log_braak_t)
pro_50_log_braak_t$braak <- sample_braak$Braak.NFT.stage[match(pro_50_log_braak_t$id, sample_braak$id)]
pro_50_log_braak_t$braak <- factor(pro_50_log_braak_t$braak, levels = 0:6)
pro_50_log_braak_t$bank <- sample_braak$bank[match(pro_50_log_braak_t$id, sample_braak$id)]
pro_50_log_braak_t$braak_type <- ifelse(pro_50_log_braak_t$braak==0,'B0',
                                        ifelse(pro_50_log_braak_t$braak%in%c(1,2),'B1',
                                               ifelse(pro_50_log_braak_t$braak%in%c(3,4),'B2','B3')))
pro_50_log_braak_t$braak_type <- factor(pro_50_log_braak_t$braak_type, levels = c('B0','B1','B2','B3'))

# abeta
pro_50_log_abeta <- pro_50_log[,colnames(pro_50_log)%in%sample_abeta$id]
pro_50_log_abeta_t <- data.frame(t(pro_50_log_abeta))
pro_50_log_abeta_t$id <- rownames(pro_50_log_abeta_t)
pro_50_log_abeta_t$abeta <- sample_abeta$A_beta_0_3[match(pro_50_log_abeta_t$id, sample_abeta$id)]
pro_50_log_abeta_t$bank <- sample_abeta$bank[match(pro_50_log_abeta_t$id, sample_abeta$id)]
pro_50_log_abeta_t$abeta_type <- ifelse(pro_50_log_abeta_t$abeta==0,'A0',
                                        ifelse(pro_50_log_abeta_t$abeta==1,'A1',
                                               ifelse(pro_50_log_abeta_t$abeta==2,'A2','A3')))
pro_50_log_abeta_t$abeta_type <- factor(pro_50_log_abeta_t$abeta_type, levels = c('A0','A1','A2','A3'))

# top信号可视化
top_gene <- c('SMOC1','SFRP1','SNRPA','SNRNP70')
top_pro <- pro_50_raw$protein[match(top_gene, pro_50_raw$gene_name)]

# adnc
for (pro in top_pro) {
  data <- pro_50_log_adnc_t[,c('id','adnc','adnc_type','bank',pro)]
  gene <- pro_50_raw$gene_name[pro_50_raw$protein==pro]
  
  y_min <- min(data[, pro], na.rm = TRUE)
  y_max <- max(data[, pro], na.rm = TRUE)
  
  my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(unique(data$adnc_type))) 
  
  p <- ggplot(data, aes(x = adnc_type, y = .data[[pro]], fill = adnc_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_bw() +
    theme(
      panel.grid = element_line(color = 'white'),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16, family = "mono"),
      plot.title = element_text(size = 16),
      legend.position = "none"
    ) +  
    labs(title = paste0(gene, '_adnc'), x = "", y = "")
  print(p)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_adnc_all.pdf'), plot = p, width = 5, height = 5)
  
  p_1 <- ggplot(data[data$bank %in% c('zju','csu'),], aes(x = adnc_type, y = .data[[pro]], fill = adnc_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title =paste0(gene,'_adnc_zju_csu'), x = "", y = "")
  print(p_1)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_adnc_zc.pdf'), plot = p_1, width = 5, height = 5)
  
  p_2 <- ggplot(data[data$bank=='pumc',], aes(x = adnc_type, y = .data[[pro]], fill = adnc_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene,'_adnc_pumc'), x = "", y = "")
  print(p_2)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_adnc_pumc.pdf'), plot = p_2, width = 5, height = 5)
}

# braak
for (i in top_pro) {
  data <- pro_50_log_braak_t[,c('id','braak','braak_type','bank',i)]
  gene <- pro_50_raw$gene_name[pro_50_raw$protein==i]
  
  y_min <- min(data[, i], na.rm = TRUE)
  y_max <- max(data[, i], na.rm = TRUE)
  
  my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(unique(data$braak_type)))
  
  p_3 <- ggplot(data, aes(x = braak_type, y = .data[[i]], fill = braak_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene, '_braak'), x = "", y = "")
  print(p_3)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_braak_all.pdf'), plot = p_3, width = 5, height = 5)
  
  
  p_4 <- ggplot(data[data$bank %in% c('zju','csu'),], aes(x = braak_type, y = .data[[i]], fill = braak_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title =paste0(gene,'_braak_zju_csu'), x = "", y = "")
  print(p_4)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_braak_zc.pdf'), plot = p_4, width = 5, height = 5)
  
  p_5 <- ggplot(data[data$bank=='pumc',], aes(x = braak_type, y = .data[[i]], fill = braak_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene,'_braak_pumc'), x = "", y = "")
  print(p_5)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_braak_pumc.pdf'), plot = p_5, width = 5, height = 5)
}

# abeta
for (i in top_pro) {
  data <- pro_50_log_abeta_t[,c('id','abeta','abeta_type','bank',i)]
  gene <- pro_50_raw$gene_name[pro_50_raw$protein==i]
  
  y_min <- min(data[, i], na.rm = TRUE)
  y_max <- max(data[, i], na.rm = TRUE)
  
  my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(unique(data$abeta_type)))
  
  p_3 <- ggplot(data, aes(x = abeta_type, y = .data[[i]], fill = abeta_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene, '_abeta'), x = "", y = "")
  print(p_3)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_abeta_all.pdf'), plot = p_3, width = 5, height = 5)
  
  
  p_4 <- ggplot(data[data$bank %in% c('zju','csu'),], aes(x = abeta_type, y = .data[[i]], fill = abeta_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title =paste0(gene,'_abeta_zju_csu'), x = "", y = "")
  print(p_4)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_abeta_zc.pdf'), plot = p_4, width = 5, height = 5)
  
  p_5 <- ggplot(data[data$bank=='pumc',], aes(x = abeta_type, y = .data[[i]], fill = abeta_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16, family = "mono"),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene,'_abeta_pumc'), x = "", y = "")
  print(p_5)
  ggsave(paste0(out_dir,'plot/box_plot/',gene,'_abeta_pumc.pdf'), plot = p_5, width = 5, height = 5)
}
