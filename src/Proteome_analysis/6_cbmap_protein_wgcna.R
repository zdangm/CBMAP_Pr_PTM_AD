# WGCNA
library(data.table)
library(impute)
library(WGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(writexl)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(limma)
library(patchwork)
library(ggprism)
library(tidyverse)
library(gground)
library(ESEA)
library(openxlsx)
library(ggtext)
library(mygene)
library(biomaRt)
library(igraph)
library(visNetwork)
library(ggraph)
library(tidygraph)
library(readr)
library(dynamicTreeCut)
library(RColorBrewer)


source('/share/home/sunly/Rscript/cbmap_pro/cbmap_586/2_sample_process.R')
out_dir <- '/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/'

##################### input data 
rownames(pro_50_raw) <- pro_50_raw$protein
pro_raw_50_im <- pro_50_raw
pro_raw_50_im <- data.frame(t(pro_raw_50_im[,colnames(pro_raw_50_im)%in%sample_adnc$id]))


##################### wgcna 
wgcna <- function(data){
  print(dim(data))
  
  gsg <- goodSamplesGenes(data)
  if (!gsg$allOK) {
    data <- data[gsg$goodSamples, gsg$goodGenes]
  }
  
  powers = c(1:20)
  sft <- pickSoftThreshold(data, 
                           networkType = "signed",
                           corFnc = "bicor",
                           powerVector = powers,
                           verbose = 5)
  # soft threshold select
  pdf(paste0(out_dir, 'plot/module/soft_threshold.pdf'), width = 8, height = 6) 
  par(mfrow = c(1,2)) 
  cex1 = 0.9
  plot(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit, signed R²",
    type = "n",
    main = "Scale Independence"
  )
  text(
    sft$fitIndices[, 1],
    -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  abline(h = 0.8, col = "blue", lty = 2)
  
  plot(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = "Mean Connectivity"
  )
  text(
    sft$fitIndices[, 1],
    sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
  )
  abline(h = 500, col = "blue", lty = 2)
  
  dev.off()
  
  sft_data <- data.frame(Power = sft$fitIndices[, 1], R2 = -sign(sft$fitIndices[,3])*sft$fitIndices[,2], mean_connectivity = sft$fitIndices[, 5])
  power_selected <- sft_data$Power[which(sft_data$R2 > 0.9 & sft_data$mean_connectivity < 500)[1]]
  print(power_selected)
  
  # 构建网络和模块识别
  cor <- WGCNA::cor
  
  net <- blockwiseModules(
    data,
    power = power_selected,
    corType = "bicor",         
    networkType = "signed",   
    TOMType = "signed",         
    deepSplit = 4,           
    pamRespectsDendro = FALSE, 
    pamStage = FALSE,
    minModuleSize = 30,         
    mergeCutHeight = 0.15, 
    numericLabels = TRUE,     
    verbose = 3
  )
  
  table(net$colors)
  
  # calculate me and kme
  MEs <- moduleEigengenes(data, net$colors)$eigengenes
  kME <- signedKME(data, MEs, corFnc = "bicor")
  
  # module distribution
  net$new_colors <- net$colors
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  
  pdf(paste0(out_dir, 'plot/module/module_tree.pdf'), width = 12, height = 8) 
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Modulecolors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  pro_num_module <- data.frame(table(moduleColors))
  
  # module color
  color_table <- data.frame(
    ModuleNumber = sort(unique(net$colors)),
    ModuleColor = labels2colors(sort(unique(net$colors)))
  )
  color_table$pro_num <- pro_num_module$Freq[match(color_table$ModuleColor, pro_num_module$moduleColors)]
  color_table$module <- paste0('M',color_table$ModuleNumber,color_table$ModuleColor)
  
  p <- ggplot(color_table, aes(x = pro_num, y = reorder(module, pro_num), fill = ModuleColor)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = pro_num), 
              hjust = -0.1,  
              size = 4) +  
    scale_fill_identity() +  
    labs(
      x = NULL,
      y = NULL,
      title = NULL
    ) +
    theme_minimal() +
    theme_bw() +
    theme(
      panel.grid=element_line(color='white'),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      plot.title = element_text(size = 14, face = "bold")
    ) +
    xlim(0, max(color_table$pro_num) * 1.1)  # 给x轴多留一点空间用于显示标签
  print(p)
  ggsave(paste0(out_dir,'plot/module/module_pro_num.pdf'), p, width = 8, height = 6, units = "in", dpi = 300)
  
  # Intramodular Connectivity
  kWithin <- intramodularConnectivity(adjacency(data, power = power_selected), net$colors)
  
  # hub protein 
  hub_pro <- chooseTopHubInEachModule(
    data,
    colorh = moduleColors,
    power = power_selected,
    type = "signed"
  )
  hub_gene <- pro_50_raw$gene_name[match(hub_pro, pro_50_raw$protein)]
  hub_gene_modules <- moduleColors[hub_pro]
  
  hub_gene_module_info <- data.frame(
    Hub_Gene = hub_gene,
    Hub_Gene_Index = hub_pro
  )
  hub_gene_module_info$Module_Color <- rownames(hub_gene_module_info)
  
  return(
    list(
      net = net,
      power =  power_selected,
      moduleColors = moduleColors,
      color_table = color_table,
      MEs = MEs,
      kME = kME,
      hub_gene = hub_gene_module_info,
      p = p
    )
  )
}

adnc_lim_all <- wgcna(pro_raw_50_im)
saveRDS(adnc_lim_all, file = "processed_data/tmp/adnc_limma_all_586.rds")
adnc_lim_all <- readRDS('processed_data/tmp/adnc_limma_all_586.rds')

hub_pro <- adnc_lim_all$hub_gene
color_table <- adnc_lim_all$color_table
color_table$module <- paste0('ME', color_table$ModuleNumber)
color_table$hub_pro <- hub_pro$Hub_Gene_Index[match(color_table$ModuleColor, hub_pro$Module_Color)]
color_table$hub_gene <- hub_pro$Hub_Gene[match(color_table$ModuleColor, hub_pro$Module_Color)]

module_protein <- list()
for (mod in color_table$ModuleNumber) {
  mod <- as.character(mod)
  protein_name <- names(adnc_lim_all$net$colors)[adnc_lim_all$net$colors==mod]
  module_protein[[mod]] <- data.frame(protein = protein_name)
  module_protein[[mod]]$module_num <- mod
  module_protein[[mod]]$module_color <- color_table$ModuleColor[match(module_protein[[mod]]$module_num, color_table$ModuleNumber)]
  module_protein[[mod]]$genename <- pro_50_raw$gene_name[match(module_protein[[mod]]$protein, pro_50_raw$protein)]
}
module_protein <- do.call(rbind,module_protein)
module_protein$module_num <- paste0('ME',module_protein$module_num)

protein_select <- pro_50_raw[pro_50_raw$protein%in%module_protein$protein[module_protein$module_num=='ME4'],]
protein_select <- protein_select[,c('gene_protein','gene_name','protein', intersect(colnames(protein_select),sample_adnc$id))]
write.table(protein_select, 'processed_data/to-ljy/module_4_protein.txt', quote = F, row.names = F, sep = '\t')
write.xlsx(module_protein, paste0(out_dir, 'table/module/module_protein.xlsx'))
write.table(color_table, paste0(out_dir, 'table/module/module_info.txt'), sep = '\t', quote = F, row.names = F)


######################### module enrichment analysis
enrich <- function(data, result, module){
  
  module_pro <- colnames(data)[result$moduleColors == module]
  module_genes <- pro_50_raw$gene_name[match(module_pro, pro_50_raw$protein)]
  gene_ids <- mapIds(org.Hs.eg.db,
                     keys = module_genes,
                     column = "ENTREZID",
                     keytype = "SYMBOL")
  
  # GO
  go_enrich <- enrichGO(
    gene = na.omit(gene_ids),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  go_enrich_BP <- enrichGO(
    gene = na.omit(gene_ids),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  go_enrich_CC <- enrichGO(
    gene = na.omit(gene_ids),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  go_enrich_MF <- enrichGO(
    gene = na.omit(gene_ids),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  ) 
  
  # KEGG
  kegg_enrich <- enrichKEGG(
    gene = na.omit(gene_ids),
    organism = "hsa",  # 人类用hsa，小鼠用mmu
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  go_result <- go_enrich@result
  if (!is.null(go_result) && nrow(go_result) > 0) {
    go_result$Gene_Symbols <- sapply(strsplit(go_result$geneID, "/"), function(x) {
      symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
      paste(symbols, collapse = ", ")
    })
  }
  
  kegg_result <- kegg_enrich@result
  if (!is.null(kegg_result) && nrow(kegg_result) > 0) {
    kegg_result$Gene_Symbols <- sapply(strsplit(kegg_result$geneID, "/"), function(x) {
      symbols <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "SYMBOL")
      paste(symbols, collapse = ", ")
    })
  }
  
  return(list(
    module_pro = module_pro,
    module_genes = module_genes,
    go_enrich = go_enrich,
    go_enrich_bp = go_enrich_BP,
    go_enrich_cc = go_enrich_CC,
    go_enrich_mf = go_enrich_MF,
    go_enrich_table = go_result,
    kegg_enrich = kegg_enrich,
    kegg_enrich_table = kegg_result))
  
}

# plot adjustment
adjust_plot_theme <- function(p) {
  p +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 80)) +
    theme(
      plot.title = element_text(size = 20),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16)
    )
}

# enrichment result organization and visualization
module_enrich_all <- list()
module_enrich_all_table <- list()
p_list <- list()

for (module in unique(adnc_lim_all$moduleColors)) {
  if (module == "grey") next 
  
  tryCatch({
    
    module_enrich_all[[module]] <- enrich(
      data = pro_raw_50_im,
      result = adnc_lim_all,
      module = module
    )
    
    # Check GO enrichment results
    if (!is.null(module_enrich_all[[module]]$go_enrich) && nrow(module_enrich_all[[module]]$go_enrich) > 0) {
      pgo <- dotplot(module_enrich_all[[module]]$go_enrich, 
                     showCategory=50,
                     title = paste0('M', adnc_lim_all$color_table$ModuleNumber[adnc_lim_all$color_table$ModuleColor==module], module, '_GO_enrichment'))
      print(pgo)
      p_list[[module]][['GO']] <- pgo
      
      pgo_bp <- dotplot(module_enrich_all[[module]]$go_enrich_bp, 
                        showCategory=50,
                        title = paste0('M', adnc_lim_all$color_table$ModuleNumber[adnc_lim_all$color_table$ModuleColor==module], module, '_GO_enrichment_BP'))
      print(pgo_bp)
      p_list[[module]][['GO_BP']] <- pgo_bp
      
      pgo_cc <- dotplot(module_enrich_all[[module]]$go_enrich_cc, 
                        showCategory=50,
                        title = paste0('M', adnc_lim_all$color_table$ModuleNumber[adnc_lim_all$color_table$ModuleColor==module], module, '_GO_enrichment_CC'))
      print(pgo_cc)
      p_list[[module]][['GO_CC']] <- pgo_cc
      
      pgo_mf <- dotplot(module_enrich_all[[module]]$go_enrich_mf, 
                        showCategory=50,
                        title = paste0('M', adnc_lim_all$color_table$ModuleNumber[adnc_lim_all$color_table$ModuleColor==module], module, '_GO_enrichment_MF'))
      print(pgo_mf)
      p_list[[module]][['GO_MF']] <- pgo_mf
      
    } else {
      message(paste("module", module, "without significant GO enrichment results"))
    }
    
    # Check KEGG enrichment results
    if (!is.null(module_enrich_all[[module]]$kegg_enrich) && nrow(module_enrich_all[[module]]$kegg_enrich) > 0) {
      pkegg <- dotplot(module_enrich_all[[module]]$kegg_enrich, 
                       showCategory=30,
                       title = paste0('M', adnc_lim_all$color_table$ModuleNumber[adnc_lim_all$color_table$ModuleColor==module], module, '_KEGG_enrichment'))
      print(pkegg)
      p_list[[module]][['KEGG']] <- pkegg
      
    } else {
      message(paste("module", module, "significant GO enrichment results"))
    }
    
  }, error = function(e) {
    message(paste("module", module, "error: ", e$message))
  })
}

# module enrichment result
module_enrich_all_table <- list()
for (name in names(module_enrich_all)) {
  module_enrich <- module_enrich_all[[name]]
  kegg <- module_enrich$kegg_enrich_table
  go <- module_enrich$go_enrich_table
  num <- adnc_lim_all$color_table$ModuleNumber[match(name, adnc_lim_all$color_table$ModuleColor)]
  module_enrich_all_table[[paste0('ME',num,'_KEGG')]] <- kegg
  module_enrich_all_table[[paste0('ME',num,'_GO')]] <- go
}

write_xlsx(module_enrich_all_table, paste0(out_dir, 'table/module/module_enrich_all.xlsx'))

# saveRDS(p_list, file = "processed_data/tmp/p_list_586.rds")
# saveRDS(module_enrich_all, file = "processed_data/tmp/wgcna_module_enrich_586.rds")
p_list <- readRDS('processed_data/tmp/p_list_586.rds')
module_enrich_all <- readRDS('processed_data/tmp/wgcna_module_enrich_586.rds')

for(module in unique(adnc_lim_all$moduleColors)) {
  ggsave(paste0('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/plot/enrichment/',module,'_GO_BP.png'),
         adjust_plot_theme(p_list[[module]][['GO_BP']]),
         width = 18,
         height = 15,
         dpi = 300)
  ggsave(paste0('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/plot/enrichment/',module,'_GO_MF.png'),
         adjust_plot_theme(p_list[[module]][['GO_MF']]),
         width = 18,
         height = 15,
         dpi = 300)
  ggsave(paste0('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/plot/enrichment/',module,'_GO_CC.png'),
         adjust_plot_theme(p_list[[module]][['GO_CC']]),
         width = 18,
         height = 15,
         dpi = 300)
  ggsave(paste0('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/plot/enrichment/',module,'_KEGG.png'),
         adjust_plot_theme(p_list[[module]][['KEGG']]),
         width = 18,
         height = 15,
         dpi = 300)
  ggsave(paste0('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/plot/enrichment/',module,'_GO.png'),
         adjust_plot_theme(p_list[[module]][['GO']]),
         width = 18,
         height = 15,
         dpi = 300)
}


# module protein information
ad_pro_meta_list <- read.xlsx('AD_pro_meta_brain_list.xlsx')
module_protein <- list()

for (mod in unique(adnc_lim_all$net$new_colors)) {
  mod_chr <- as.character(mod)
  
  proteins_in_mod <- names(adnc_lim_all$net$new_colors)[adnc_lim_all$net$new_colors == mod]
  
  if (length(proteins_in_mod) == 0) next
  
  module_protein[[mod_chr]] <- data.frame(
    pro = proteins_in_mod
  )
  
  module_protein[[mod_chr]]$gene <- pro_50_raw$gene_name[
    match(module_protein[[mod_chr]]$pro, pro_50_raw$protein)
  ]
  module_protein[[mod_chr]]$module <- mod_chr
}

module_protein <- do.call(rbind, module_protein)
module_protein$module <- paste0('ME', module_protein$module)

me4_reported <- intersect(module_protein$gene[module_protein$module=='ME4'], ad_pro_meta_list$gene[ad_pro_meta_list$FDR<0.05])
me4_notreported <- setdiff(module_protein$gene[module_protein$module=='ME4'], me4_reported)

me10_reported <- intersect(module_protein$gene[module_protein$module=='ME10'], ad_pro_meta_list$gene[ad_pro_meta_list$FDR<0.05])
me10_notreported <- setdiff(module_protein$gene[module_protein$module=='ME10'], me10_reported)

me15_reported <- intersect(module_protein$gene[module_protein$module=='ME15'], ad_pro_meta_list$gene[ad_pro_meta_list$FDR<0.05])
me15_notreported <- setdiff(module_protein$gene[module_protein$module=='ME15'], me15_reported)

sample_id <- colnames(pro_50_raw)[colnames(pro_50_raw)%in%sample_adnc$id]
module_protein_select <- merge(module_protein, pro_50_raw[,c('protein', sample_id)], by.x = 'pro', by.y = 'protein')
module_protein_select <- module_protein_select[!module_protein_select$module=='ME0',]


############################### module trait analysis

# module & disease
module_trait_cor <- function(data){
  me <- data[,2:ncol(data)]
  
  trait <- data.frame(
    id = sample_infor$id[match(rownames(me),sample_infor$id)],
    adnc = sample_infor$adnc_num[match(rownames(me),sample_infor$id)],
    age = sample_infor$age[match(rownames(me),sample_infor$id)],
    sex = sample_infor$sex_male[match(rownames(me),sample_infor$id)],
    braak = sample_infor$Braak.NFT.stage[match(rownames(me),sample_infor$id)],
    abeta = sample_infor$A_beta_0_3[match(rownames(me),sample_infor$id)],
    cscore = sample_infor$C_cerad_0_3[match(rownames(me),sample_infor$id)],
    bank = sample_infor$bank[match(rownames(me),sample_infor$id)],
    rin = sample_infor$RIN[match(rownames(me),sample_infor$id)],
    pmd = sample_infor$PMD[match(rownames(me),sample_infor$id)],
    part = sample_infor$part[match(rownames(me),sample_infor$id)],
    late = sample_infor$late[match(rownames(me),sample_infor$id)],
    artag = sample_infor$artag[match(rownames(me),sample_infor$id)],
    lbd = sample_infor$lbd[match(rownames(me),sample_infor$id)],
    cvd = sample_infor$cvd[match(rownames(me),sample_infor$id)]
  )
  
  # adnc
  me <- t(me)
  trait_adnc <- trait[!is.na(trait$adnc),]
  me_adnc <- me[,colnames(me)%in%trait_adnc$id]
  me_adnc <- me_adnc[,match(trait_adnc$id,colnames(me_adnc))]
  design <- model.matrix( ~ adnc + age + sex + pmd + rin, data = trait_adnc)
  fit <- lmFit(me_adnc, design)
  fit <- eBayes(fit)
  adnc_num_lim_re <- topTable(fit, coef = "adnc", number = Inf, adjust.method = "BH")
  adnc_num_lim_re$me <- rownames(adnc_num_lim_re)
  adnc_num_lim_re$trait <- 'ADNC'
  
  # age
  age_lim_re <- topTable(fit, coef = "age", number = Inf, adjust.method = "BH")
  age_lim_re$me <- rownames(age_lim_re)
  age_lim_re$trait <- 'Age'
  
  # braak
  trait_braak <- trait[trait$id%in%sample_braak$id,]
  me_braak <- me[,colnames(me)%in%trait_braak$id]
  me_braak <- me_braak[,match(trait_braak$id,colnames(me_braak))]
  design <- model.matrix( ~ braak + age + sex + pmd + rin, data = trait_braak)
  fit <- lmFit(me_braak, design)
  fit <- eBayes(fit)
  braak_num_lim_re <- topTable(fit, coef = "braak", number = Inf, adjust.method = "BH")
  braak_num_lim_re$me <- rownames(braak_num_lim_re)
  braak_num_lim_re$trait <- 'Braak'
  
  # abeta
  trait_abeta <- trait[trait$id%in%sample_abeta$id,]
  me_abeta <- me[,colnames(me)%in%trait_abeta$id]
  me_abeta <- me_abeta[,match(trait_abeta$id,colnames(me_abeta))]
  design <- model.matrix( ~ abeta + age + sex + pmd + rin, data = trait_abeta)
  fit <- lmFit(me_abeta, design)
  fit <- eBayes(fit)
  abeta_num_lim_re <- topTable(fit, coef = "abeta", number = Inf, adjust.method = "BH")
  abeta_num_lim_re$me <- rownames(abeta_num_lim_re)
  abeta_num_lim_re$trait <- 'Abeta'
  
  # Cerad
  trait_cscore <- trait[trait$id%in%sample_cscore$id,]
  me_cscore <- me[,colnames(me)%in%trait_cscore$id]
  me_cscore <- me_cscore[,match(trait_cscore$id,colnames(me_cscore))]
  design <- model.matrix( ~ cscore + age + sex + pmd + rin, data = trait_cscore)
  fit <- lmFit(me_cscore, design)
  fit <- eBayes(fit)
  cscore_num_lim_re <- topTable(fit, coef = "cscore", number = Inf, adjust.method = "BH")
  cscore_num_lim_re$me <- rownames(cscore_num_lim_re)
  cscore_num_lim_re$trait <- 'Cerad'
  
  # part
  trait_part <- trait[trait$id%in%sample_part$id,]
  me_part <- me[,colnames(me)%in%trait_part$id]
  me_part <- me_part[,match(trait_part$id,colnames(me_part))]
  design <- model.matrix( ~ part + age + sex + pmd + rin, data = trait_part)
  fit <- lmFit(me_part, design)
  fit <- eBayes(fit)
  part_num_lim_re <- topTable(fit, coef = "part", number = Inf, adjust.method = "BH")
  part_num_lim_re$me <- rownames(part_num_lim_re)
  part_num_lim_re$trait <- 'PART'
  
  # late
  trait_late <- trait[trait$id%in%sample_late$id,]
  me_late <- me[,colnames(me)%in%trait_late$id]
  me_late <- me_late[,match(trait_late$id,colnames(me_late))]
  design <- model.matrix( ~ late + age + sex + pmd + rin, data = trait_late)
  fit <- lmFit(me_late, design)
  fit <- eBayes(fit)
  late_num_lim_re <- topTable(fit, coef = "late", number = Inf, adjust.method = "BH")
  late_num_lim_re$me <- rownames(late_num_lim_re)
  late_num_lim_re$trait <- 'LATE'
  
  # artag
  trait_artag <- trait[trait$id%in%sample_artag$id,]
  me_artag <- me[,colnames(me)%in%trait_artag$id]
  me_artag <- me_artag[,match(trait_artag$id,colnames(me_artag))]
  design <- model.matrix( ~ artag + age + sex + pmd + rin, data = trait_artag)
  fit <- lmFit(me_artag, design)
  fit <- eBayes(fit)
  artag_num_lim_re <- topTable(fit, coef = "artag", number = Inf, adjust.method = "BH")
  artag_num_lim_re$me <- rownames(artag_num_lim_re)
  artag_num_lim_re$trait <- 'ARTAG'
  
  # lbd
  trait_lbd <- trait[trait$id%in%sample_lbd$id,]
  me_lbd <- me[,colnames(me)%in%trait_lbd$id]
  me_lbd <- me_lbd[,match(trait_lbd$id,colnames(me_lbd))]
  design <- model.matrix( ~ lbd + age + sex + pmd + rin, data = trait_lbd)
  fit <- lmFit(me_lbd, design)
  fit <- eBayes(fit)
  lbd_num_lim_re <- topTable(fit, coef = "lbd", number = Inf, adjust.method = "BH")
  lbd_num_lim_re$me <- rownames(lbd_num_lim_re)
  lbd_num_lim_re$trait <- 'LBD'
  
  # cvd
  trait_cvd <- trait[trait$id%in%sample_cvd$id,]
  me_cvd <- me[,colnames(me)%in%trait_cvd$id]
  me_cvd <- me_cvd[,match(trait_cvd$id,colnames(me_cvd))]
  design <- model.matrix( ~ cvd + age + sex + pmd + rin, data = trait_cvd)
  fit <- lmFit(me_cvd, design)
  fit <- eBayes(fit)
  cvd_num_lim_re <- topTable(fit, coef = "cvd", number = Inf, adjust.method = "BH")
  cvd_num_lim_re$me <- rownames(cvd_num_lim_re)
  cvd_num_lim_re$trait <- 'CVD'
  
  trait_table <- rbind(adnc_num_lim_re, age_lim_re, braak_num_lim_re, abeta_num_lim_re, cscore_num_lim_re, part_num_lim_re, late_num_lim_re, artag_num_lim_re, lbd_num_lim_re, cvd_num_lim_re)
  
  return(
    list(
      adnc_num_lim_re = adnc_num_lim_re,
      age_lim_re = age_lim_re,
      braak_num_lim_re = braak_num_lim_re,
      trait_table = trait_table
    )
  )
}

module_trait_all <- module_trait_cor(adnc_lim_all$MEs)

module_trait_all$trait_table$pmarker <- ifelse(module_trait_all$trait_table$adj.P.Val < 0.05, "*", '')
module_trait_all$trait_table$me <- factor(module_trait_all$trait_table$me, levels = paste0('ME', 1:length(unique(module_trait_all$trait_table$me))))
module_trait_all$trait_table$trait <- factor(module_trait_all$trait_table$trait, levels = c('Age','LATE','ARTAG','PART','CVD','LBD','Cerad','Braak','Abeta','ADNC'))

p1 <- ggplot(module_trait_all$trait_table, aes(x = me, y = trait, fill = -log10(P.Value))) +
  geom_tile(color = "white") +
  geom_text(aes(label = pmarker), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "#fd8d3c", name = "-log10(P)") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA), 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 18),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = NULL, y = NULL, title = NULL)
p1

write.xlsx(module_trait_all$trait_table, paste0(out_dir, 'table/module/module_trait_asso.xlsx'))


# module 细胞类型特异性分析
cell_type_marker <- read.xlsx('cell_type_markers.xlsx')
colnames(cell_type_marker) <- c('Astro','Micro','Neuro','Oligo','Endo')

module_marker_enrichment <- function(data){
  
  all_proteins <- data.frame(pro = names(data$net$new_colors))
  all_proteins$gene <- pro_raw$gene_name[match(all_proteins$pro, pro_raw$protein)]
  
  module_colors <- unique(data$net$new_colors)
  module_enrichment_results <- list()
  for (module in module_colors) {
    if (module == 0) next
    
    module_proteins <- names(data$net$new_colors)[data$net$new_colors == module]
    module_proteins <- all_proteins$gene[all_proteins$pro %in% module_proteins]
    
    for (cell_type in colnames(cell_type_marker)) {
      markers <- cell_type_marker[[cell_type]]
      markers <- markers[!is.na(markers)]
      
      # a: in module and marker genes
      a <- sum(module_proteins %in% markers)
      
      # b: in module but not marker genes
      b <- sum(!(module_proteins %in% markers))
      
      # c: not in module but are marker genes
      other_proteins <- setdiff(all_proteins$gene, module_proteins)
      c <- sum(other_proteins %in% markers)
      
      # d: not in module and not marker genes
      d <- sum(!(other_proteins %in% markers))
      
      # Fisher test
      mat <- matrix(c(a, b, c, d), nrow = 2)
      ft <- fisher.test(mat)
      
      # result save
      module_enrichment_results[[paste0("M", module,data$color_table$ModuleColor[data$color_table$ModuleNumber==module], "_", cell_type)]] <- list(
        module = module,
        cell_type = cell_type,
        p.value = ft$p.value,
        odds_ratio = ft$estimate,
        overlap = a,
        module_size = length(module_proteins),
        marker_size = length(markers)
      )
    }
  }
  
  enrich_df <- do.call(rbind, lapply(module_enrichment_results, as.data.frame))
  enrich_df <- data.frame(enrich_df)
  enrich_df <- enrich_df[order(enrich_df$module),]
  enrich_df$FDR <- p.adjust(enrich_df$p.value, method = "BH")
  enrich_df$color <- data$color_table$ModuleColor[match(enrich_df$module,data$color_table$ModuleNumber)]
  enrich_df$module_color <- paste0('M',enrich_df$module,enrich_df$color)
  enrich_df$module_color <- factor(enrich_df$module_color, levels = data$color_table$module[-1])
  enrich_df$module_num <- paste0('ME', enrich_df$module)
  
  # plot data preparation
  plot_df <- enrich_df %>%
    mutate(
      logP = -log10(p.value),  # 避免 -Inf
      module_num = factor(module_num, levels = paste0('ME',1:length(unique(enrich_df$module_num)))),
      cell_type = factor(cell_type)
    )
  plot_df$signif <- ifelse(plot_df$FDR < 0.05, "*", "")
  
  p <- ggplot(plot_df, aes(x = module_num, y = cell_type, fill = logP)) +
    geom_tile(color = "white") +
    geom_text(aes(label = signif), color = "black", size = 5) +
    scale_fill_gradient(low = "white", high = "#fd8d3c", name = "-log10(P)") +
    theme_minimal(base_size = 12) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 18),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, margin = margin(b = 4))
    )
  print(p)
  
  return(list(result = enrich_df, plot = p))
}

cell_maker_all <- module_marker_enrichment(adnc_lim_all)
p2 <- cell_maker_all$plot

write.xlsx(cell_maker_all$result, paste0(out_dir,'table/module/cell_typ.xlsx'))

# module and NFT/Abeta plaque protein
abeta_tau_marker <- read.xlsx('abeta_tau.xlsx')
colnames(abeta_tau_marker) <- c('NFT','pTau','pTau NFT','Abeta plaque','Abeta con-plaque')

abeta_tau_marker_enrichment <- function(data){
  
  all_proteins <- data.frame(pro = names(data$net$new_colors))
  all_proteins$gene <- pro_raw$gene_name[match(all_proteins$pro, pro_raw$protein)]
  
  module_colors <- unique(data$net$new_colors)
  module_enrichment_results <- list()
  for (module in module_colors) {
    if (module == 0) next
    
    module_proteins <- names(data$net$new_colors)[data$net$new_colors == module]
    module_proteins <- all_proteins$gene[all_proteins$pro %in% module_proteins]
    
    for (type in colnames(abeta_tau_marker)) {
      markers <- abeta_tau_marker[[type]]
      markers <- markers[!is.na(markers)]
      
      a <- sum(module_proteins %in% markers)
      
      b <- sum(!(module_proteins %in% markers))
      
      other_proteins <- setdiff(all_proteins$gene, module_proteins)
      c <- sum(other_proteins %in% markers)
      
      d <- sum(!(other_proteins %in% markers))
      
      mat <- matrix(c(a, b, c, d), nrow = 2)
      ft <- fisher.test(mat)
      
      # result save 
      module_enrichment_results[[paste0("M", module,data$color_table$ModuleColor[data$color_table$ModuleNumber==module], "_", type)]] <- list(
        module = module,
        type = type,
        p.value = ft$p.value,
        odds_ratio = ft$estimate,
        overlap = a,
        module_size = length(module_proteins),
        marker_size = length(markers)
      )
    }
  }
  
  enrich_df <- do.call(rbind, lapply(module_enrichment_results, as.data.frame))
  enrich_df <- data.frame(enrich_df)
  enrich_df$FDR <- p.adjust(enrich_df$p.value, method = "BH")
  enrich_df$color <- data$color_table$ModuleColor[match(enrich_df$module,data$color_table$ModuleNumber)]
  enrich_df$module_color <- paste0('M',enrich_df$module,enrich_df$color)
  enrich_df$module_color <- factor(enrich_df$module_color, levels = data$color_table$module[-1])
  enrich_df$module_num <- paste0('ME', enrich_df$module)
  enrich_df$type <- factor(enrich_df$type, levels = c('pTau','pTau NFT','NFT','Abeta con-plaque','Abeta plaque'))
  
  # plot data prepare
  plot_df <- enrich_df %>%
    mutate(
      logP = -log10(p.value), 
      module_num = factor(module_num, levels = paste0('ME',1:length(unique(enrich_df$module_num)))),
      cell_type = factor(type)
    )
  plot_df$signif <- ifelse(plot_df$FDR < 0.05, "*", "")
  
  base_labels <- paste0("ME", 1:length(unique(enrich_df$module_num)))
  special_positions <- c(4, 11, 16) 
  styled_labels <- ifelse(1:length(unique(enrich_df$module_num)) %in% special_positions,
                          paste0("<span style='color:red; font-weight:bold;'>", base_labels, "</span>"),
                          base_labels)
  
  p <- ggplot(plot_df, aes(x = module_num, y = type, fill = logP)) +
    geom_tile(color = "white") +
    geom_text(aes(label = signif), color = "black", size = 5) +
    scale_fill_gradient(low = "white", high = "#fd8d3c", name = "-log10(P)") +
    scale_x_discrete(labels = styled_labels) +  # Add this line
    theme_minimal(base_size = 12) +
    labs(
      title = NULL,
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.x = element_markdown(size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 18),
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, margin = margin(b = 4))
    )
  print(p)
  
  return(list(result = enrich_df, plot = p))
}

abeta_tau_maker_all <- abeta_tau_marker_enrichment(adnc_lim_all)

p3 <- abeta_tau_maker_all$plot

p_combined <- (p1 / p2 / p3) + plot_layout(heights = c(2, 1, 1))
p_combined

ggsave(paste0(out_dir, 'plot/heatmap/trait_cell_at_com_586.pdf'),
       p_combined,
       width = 8,
       height = 10)
write.xlsx(abeta_tau_maker_all$result, paste0(out_dir,'table/module/abeta_tau.xlsx'))


######################### highlight enrichment result of three(4,11,16) module
pal <- c(
  "KEGG" = '#c3e1e6',
  "MF" = '#f3dfb7',
  "CC" = '#dcc6dc',
  "BP" = '#96c38e'
)

# me4
ego_readable_4 <- setReadable(module_enrich_all$yellow$go_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ekegg_readable_4 <- setReadable(module_enrich_all$yellow$kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
GO_4 <- as.data.frame(ego_readable_4)
KEGG_4 <- as.data.frame(ekegg_readable_4)

top_bp_4 <- GO_4 %>%
  filter(ONTOLOGY == "BP") %>%
  top_n(5, wt = -pvalue)

top_cc_4 <- GO_4 %>%
  filter(ONTOLOGY == "CC") %>%
  top_n(5, wt = -pvalue)

top_mf_4 <- GO_4 %>%
  filter(ONTOLOGY == "MF") %>%
  top_n(5, wt = -pvalue)

top_kegg_4 <- KEGG_4 %>%
  top_n(11, wt = -pvalue) %>%
  mutate(ONTOLOGY = "KEGG")

use_pathway_4 <- bind_rows(top_bp_4, top_cc_4, top_mf_4, top_kegg_4) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
  arrange(ONTOLOGY, -pvalue) %>%
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
  tibble::rowid_to_column("index") %>%
  mutate(module = 'ME4')

# me11
ego_readable_11 <- setReadable(module_enrich_all$greenyellow$go_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ekegg_readable_11 <- setReadable(module_enrich_all$greenyellow$kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
GO_11 <- as.data.frame(ego_readable_11)
KEGG_11 <- as.data.frame(ekegg_readable_11)

top_bp_11 <- GO_11 %>%
  filter(ONTOLOGY == "BP") %>%
  top_n(4, wt = -pvalue)

top_cc_11 <- GO_11 %>%
  filter(ONTOLOGY == "CC") %>%
  top_n(4, wt = -pvalue)

top_mf_11 <- GO_11 %>%
  filter(ONTOLOGY == "MF") %>%
  top_n(4, wt = -pvalue)

top_kegg_11 <- KEGG_11 %>%
  top_n(2, wt = -pvalue) %>%
  mutate(ONTOLOGY = "KEGG")

use_pathway_11 <- bind_rows(top_bp_11, top_cc_11, top_mf_11, top_kegg_11) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF","KEGG")))) %>%
  arrange(ONTOLOGY, -pvalue) %>%
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = Description)) %>%  # ← 每个组内部逆序
  tibble::rowid_to_column("index") %>%
  mutate(module = 'ME11')

# me16
ego_readable_16 <- setReadable(module_enrich_all$lightcyan$go_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
ekegg_readable_16 <- setReadable(module_enrich_all$lightcyan$kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
GO_16 <- as.data.frame(ego_readable_16)
KEGG_16 <- as.data.frame(ekegg_readable_16)

top_bp_16 <- GO_16 %>%
  filter(ONTOLOGY == "BP") %>%
  top_n(5, wt = -pvalue)

top_cc_16 <- GO_16 %>%
  filter(ONTOLOGY == "CC") %>%
  top_n(5, wt = -pvalue)

top_mf_16 <- GO_16 %>%
  filter(ONTOLOGY == "MF") %>%
  top_n(4, wt = -pvalue)

top_kegg_16 <- KEGG_16 %>%
  top_n(2, wt = -pvalue) %>%
  mutate(ONTOLOGY = "KEGG")

use_pathway_16 <- bind_rows(top_bp_16, top_cc_16, top_mf_16, top_kegg_16) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF",'KEGG')))) %>%
  arrange(ONTOLOGY, -pvalue) %>%
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  tibble::rowid_to_column("index") %>%
  mutate(module = 'ME16')


# plot
xmax1 <- max(-log10(use_pathway_4$pvalue))
xmax2 <- max(-log10(use_pathway_11$pvalue))
xmax3 <- max(-log10(use_pathway_16$pvalue))

xaxis_max <- max(c(xmax1, xmax2, xmax3)) + 1

global_count_min <- min(
  min(use_pathway_4$Count),
  min(use_pathway_11$Count),
  min(use_pathway_16$Count)
)

global_count_max <- max(
  max(use_pathway_4$Count),
  max(use_pathway_11$Count),
  max(use_pathway_16$Count)
)

pmod4 <- use_pathway_4 %>% 
  ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5.5
  ) +
  geom_point(
    aes(x = -0.02 * xaxis_max, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.02 * xaxis_max, label = Count)
  ) +
  scale_size_continuous(name = "Count", range = c(3, 12), limits = c(global_count_min, global_count_max)) +
  labs(x = NULL, y = NULL, title = 'ME4') +
  scale_fill_manual(name = "Category", values = pal, drop = FALSE) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    plot.title = element_text(size = 14, margin = margin(b = 2)),
    plot.margin = margin(t = 2, r = 5, b = 0, l = 5),
    axis.text.x = element_blank(),    
    axis.ticks.x = element_blank(),  
    axis.line.x = element_blank()    
  )
pmod4

pmod11 <- use_pathway_11 %>% 
  ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.05, label = Description),
    hjust = 0, size = 5.5
  ) +
  geom_point(
    aes(x = -0.004 * xaxis_max, size = Count),
    shape = 21
  ) +
  geom_text(
    aes(x = -0.004 * xaxis_max, label = Count)
  ) +
  labs(x = NULL, y = NULL, title = 'ME11') +
  scale_size_continuous(name = "Count", range = c(5, 12), limits = c(global_count_min, global_count_max)) +
  scale_fill_manual(name = "Category", values = pal, drop = FALSE) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'none',
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(), 
    axis.line.x = element_blank(),
    plot.title = element_text(size = 14, margin = margin(b = 2)),
    plot.margin = margin(t = 2, r = 5, b = 0, l = 5),
  )
pmod11

fake_kegg <- data.frame(
  pvalue = 1,
  index = -999,  # 不影响实际图层
  Description = "Fake KEGG",
  Count = 0,
  ONTOLOGY = factor("KEGG", levels = levels(use_pathway_16$ONTOLOGY)),
  module = "ME16"
)

use_pathway_16_fake <- bind_rows(fake_kegg, use_pathway_16)
use_pathway_16_fake$index <- 1:nrow(use_pathway_16_fake)

use_pathway_16_fake <- use_pathway_16_fake %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF",'KEGG')))) %>%
  arrange(ONTOLOGY, -pvalue) %>%
  group_by(ONTOLOGY) %>%
  mutate(Description = factor(Description, levels = Description)) %>% 
  mutate(module = 'ME16')

pmod16 <- use_pathway_16_fake %>% 
  ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
  geom_round_col(
    aes(y = Description), width = 0.6, alpha = 0.8
  ) +
  geom_text(
    aes(x = 0.05, label = Description, 
        color = ifelse(Description == "Fake KEGG", NA, "black")),
    hjust = 0, size = 5.5,
    show.legend = FALSE
  ) +
  geom_point(
    aes(x = -0.017 * xaxis_max, size = Count),
    shape = 21,
    alpha = 1,
    show.legend = TRUE
  ) +
  geom_point(
    data = filter(use_pathway_16_fake, Description == "Fake KEGG"),
    aes(x = -0.017 * xaxis_max, size = Count),
    shape = 21,
    alpha = 0,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = -0.017 * xaxis_max, label = Count,
        color = ifelse(Description == "Fake KEGG", NA, "black")),
    show.legend = FALSE
  ) +
  labs(y = NULL, title = 'ME16') +
  scale_size_continuous(name = "Count", range = c(5, 12), limits = c(global_count_min, global_count_max)) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    data = data.frame(xaxis_max = xaxis_max),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(name = "Category", 
                    values = pal[c("BP", "CC", "MF", "KEGG")], 
                    breaks = c("BP", "CC", "MF", "KEGG"),
                    limits = c("BP", "CC", "MF", "KEGG"),  
                    drop = FALSE) +
  theme_prism() +
  theme(
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(size = 14, margin = margin(b = 2)),
    plot.margin = margin(t = 2, r = 5, b = 0, l = 5),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 15)
  ) +
  guides(
    size = guide_legend(override.aes = list(shape = 21)),
    fill = guide_legend(override.aes = list(shape = NA, fill = pal[c("BP", "CC", "MF", "KEGG")]))
  ) +
  scale_color_identity()

pmod16

pmod <- pmod4 + pmod11 + pmod16 + plot_layout(ncol = 1, heights = c(1.5, 1, 1))
pmod
ggsave(paste0(out_dir, 'plot/enrichment/module_enrich.pdf'),
       plot = pmod,
       width = 12,
       height = 15)


############################## PPI
ppi_brain <- read.xlsx('ppi_list.xlsx')
resa <- queryMany(ppi_brain$protein_A, scopes="uniprot", fields="symbol", species="human")
genea <- data.frame(resa[, c("query", "symbol")])
resb <- queryMany(ppi_brain$protein_B, scopes="uniprot", fields="symbol", species="human")
geneb <- data.frame(resb[, c("query", "symbol")])
ppi_brain$gene_A <- genea$symbol[match(ppi_brain$protein_A, genea$query)]
ppi_brain$gene_B <- geneb$symbol[match(ppi_brain$protein_B, geneb$query)]

ppi_string <- as.data.frame(fread('string_data/9606.protein.physical.links.v12.0.min900.onlyAB.txt'))
ppi_string$protein1 <- sub("9606.", "", ppi_string$protein1)
ppi_string$protein2 <- sub("9606.", "", ppi_string$protein2)
gene_id <- unique(c(ppi_string$protein1, ppi_string$protein2))
res <- queryMany(gene_id, scopes = "ensembl.protein", fields = "symbol", species = "human")
gene <- data.frame(res[,c("query", "symbol")])
ppi_string$gene_A <- gene$symbol[match(ppi_string$protein1, gene$query)]
ppi_string$gene_B <- gene$symbol[match(ppi_string$protein2, gene$query)]

ppi_293 <- read_tsv('string_data/BioPlex_293T_Network_10K_Dec_2019.tsv')
ppi_293 <- ppi_293[ppi_293$pInt>0.9,]
colnames(ppi_293)[5:6] <- c('gene_A','gene_B')
ppi_116 <- read_tsv('string_data/BioPlex_HCT116_Network_5.5K_Dec_2019.tsv')
ppi_116 <- ppi_116[ppi_116$pInt>0.9,]
colnames(ppi_116)[5:6] <- c('gene_A','gene_B')

ppi_all <- rbind(ppi_brain[,11:12], ppi_string[,4:5], ppi_293[,5:6], ppi_116[5:6])
ppi_all <- na.omit(ppi_all)
ppi_all <- unique(ppi_all)

write.table(ppi_all, paste0(out_dir,'table/module/ppi_list_all.txt'), sep = '\t', quote = F, row.names = F)

# ppi pengjunmin
analyze_ppi_modules_louvain <- function(data, min_edges = 5, minClusterSize = 5,
                                        pval_cutoff = 0.05, qval_cutoff = 0.2, ont = "ALL") {
  
  modules <- unique(data$net$new_colors)
  
  ppi_modules <- list()
  ppi_edges <- list()
  enrichment_results <- list()
  modularity_score <- list()
  module_sizes <- list()
  
  for (mod in modules) {
    
    message("Processing WGCNA module: ", mod)
    
    if (mod == 0) next
    
    module_proteins <- names(data$net$new_colors)[data$net$new_colors == mod]
    module_genes <- pro_50_raw$gene_name[match(module_proteins, pro_50_raw$protein)]
    mod_edges <- subset(ppi_all, gene_A %in% module_genes & gene_B %in% module_genes)
    ppi_edges[[mod]] <- mod_edges
    
    if (nrow(mod_edges) >= min_edges) {
      
      g <- igraph::graph_from_data_frame(mod_edges, directed = FALSE)
      cl <- igraph::cluster_louvain(g)
      
      modularity_score[[mod]] <- modularity(cl)
      module_sizes[[mod]] <- sizes(cl) %>% 
        as.data.frame() %>%
        setNames(c("Module", "Size")) %>%
        arrange(desc(Size))
      
      cat("Modularity:", modularity_score[[mod]], "\n")
      
      gene_names <- igraph::V(g)$name
      submods <- igraph::membership(cl)
      
      module_result <- data.frame(Gene = gene_names, SubModule = submods)
      
      all_genes <- module_genes[!is.na(module_genes)]
      unassigned_genes <- setdiff(all_genes, gene_names)
      
      if (length(unassigned_genes) > 0) {
        module_result <- rbind(module_result, data.frame(Gene = unassigned_genes, SubModule = 0))
      }
      
      ppi_modules[[mod]] <- module_result
      
      # GO 
      enrich_mod <- list()
      for (submod in unique(module_result$SubModule)) {
        
        if (submod == 0) next
        
        sub_genes <- module_result$Gene[module_result$SubModule == submod]
        
        if (length(sub_genes) < minClusterSize) next
        
        gene_entrez <- suppressMessages(
          clusterProfiler::bitr(sub_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
        )
        
        if (nrow(gene_entrez) >= 5) {
          ego <- clusterProfiler::enrichGO(gene = gene_entrez$ENTREZID,
                                           OrgDb = org.Hs.eg.db,
                                           ont = ont,
                                           pAdjustMethod = "BH",
                                           pvalueCutoff = pval_cutoff,
                                           qvalueCutoff = qval_cutoff,
                                           readable = TRUE)
          enrich_mod[[paste0("submod_", submod)]] <- ego
        }
      }
      
      enrichment_results[[mod]] <- enrich_mod
      
    } else {
      warning("Module ", mod, " has too few PPI edges and was skipped")
    }
  }
  
  return(
    list(
      ppi_modules = ppi_modules, 
      enrichment = enrichment_results, 
      ppi_edges = ppi_edges,
      modularity_score = modularity_score,
      module_sizes = module_sizes
    )
  )
  
}

adnc_ppi_result_louvain <- analyze_ppi_modules_louvain(data = adnc_lim_all)
# saveRDS(adnc_ppi_result_louvain, file = "processed_data/tmp/adnc_ppi_louvain_586.rds")
# adnc_ppi_result_louvain <- readRDS('processed_data/tmp/adnc_ppi_louvain.rds')

submod_enrich <- list()
for (i in c(4,11,16)) {
  module_enrich <- adnc_ppi_result_louvain$enrichment[[i]]
  list_name <- names(module_enrich)
  for(name in list_name){
    table <- module_enrich[[name]]@result
    submod_enrich[[paste0('ME',i,'_',name)]] <- table
  }
}
write.xlsx(submod_enrich,paste0(out_dir, 'table/module/submod_enrich.xlsx'))

# ppi_可视化
adnc_num_lim_re <- as.data.frame(fread(paste0(out_dir,'table/dea_result/cbmap_adnc_lim_re.txt')))

select_module <- c(4,11,16)

ppi_mod_edge <- list()
ppi_mod_node <- list()

for (mod in select_module) {
  
  data <- adnc_ppi_result_louvain$ppi_modules[[mod]]
  module_size <- data.frame(table(data$SubModule))
  colnames(module_size) <- c('SubModule','size')
  data$module_size <- module_size$size[match(data$SubModule, module_size$SubModule)]
  data <- data[data$SubModule != 0 & data$module_size >= 5, ]
  data$type <- ifelse(data$Gene %in% adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val < 0.05], 'DEP', 'not-DEP')
  data$type <- as.factor(data$type)
  
  if (!is.null(data) && nrow(data) > 0) {
    data$SubModule <- paste0('m', data$SubModule)
    ppi_mod_node[[mod]] <- data
    ppi_mod_node[[mod]]$wgcna_module <- paste0('ME', mod)
    ppi_mod_node[[mod]]$ppi_module <- paste0(ppi_mod_node[[mod]]$wgcna_module, '_', ppi_mod_node[[mod]]$SubModule)
  }
  
  if (is.null(data) || nrow(data) == 0) next
  
  # 1. 显式地把 SubModule 转成字符，再转因子，避免乱序
  ordered_levels <- paste0("m", sort(as.numeric(gsub("m", "", unique(data$SubModule)))))
  data$SubModule <- factor(data$SubModule, levels = ordered_levels)
  submodule_levels <- levels(data$SubModule)
  n_submodules <- length(submodule_levels)
  
  # 2. 构建具名颜色向量，确保 SubModule 与颜色一一对应
  if (n_submodules > 12) {
    palette_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_submodules)
  } else {
    palette_colors <- brewer.pal(max(n_submodules, 3), "Set3")[1:n_submodules]
  }
  
  names(palette_colors) <- submodule_levels  # 关键：确保名字与 levels 一致
  
  # 3. 构建图
  mod_edges <- subset(ppi_all, gene_A %in% data$Gene & gene_B %in% data$Gene)
  ppi_mod_edge[[mod]] <- mod_edges
  
  ppi_graph <- tbl_graph(
    nodes = data, 
    edges = mod_edges %>% rename(from = gene_A, to = gene_B),
    directed = FALSE
  )
  
  print(head(as_tibble(ppi_graph, active = "nodes")))
  # 4. 绘图
  p <- ggraph(ppi_graph, layout = "fr") + 
    geom_edge_link(color = "grey60") +
    
    geom_node_point(
      aes(fill = SubModule, shape = type), 
      size = 4,
      stroke = 0.5,
      color = "grey70" # 不遮盖 fill
    ) +
    
    geom_node_text(
      aes(label = Gene),  
      size = 3,
      repel = TRUE,
      family = "Arial",
      box.padding = unit(0.4, "lines"),
      point.padding = unit(0.3, "lines")
    ) +
    
    scale_shape_manual(
      values = c("DEP" = 24, "not-DEP" = 21),
      name = "Protein type"
    ) +
    
    scale_fill_manual(
      values = palette_colors,
      name = "SubModule"
    )  +
    
    theme_void() +
    
    theme(
      legend.position = 'right',
      legend.box = "vertical",
      legend.spacing.y = unit(0.5, "cm"),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    
    labs(
      title = paste0('PPI of ', adnc_lim_all$color_table$module[adnc_lim_all$color_table$ModuleNumber == mod]),
      fill = "SubModule"  # 显式指定 fill 图例名
    ) +
    
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 5)),
      shape = guide_legend(override.aes = list(fill = "grey70", size = 5))
    )
  
  print(p)
  
  ggsave(paste0(out_dir,'plot/ppi/m_',mod,'.png'),
         p,
         width = 15,
         height = 10)
}

ppi_mod_edge_all <- do.call(rbind, ppi_mod_edge)
ppi_mod_edge_all <- ppi_mod_edge_all %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(gene_A, gene_B)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)
ppi_mod_node_all <- do.call(rbind, ppi_mod_node)
ppi_mod_node_all$hub_gene <- ifelse(ppi_mod_node_all$Gene%in%adnc_lim_all$hub_gene$Hub_Gene,'Yes','No')

# 绘制三个module的ppi信息于一张图上
ppi_graph <- tbl_graph(
  nodes = ppi_mod_node_all, 
  edges = ppi_mod_edge_all %>% rename(from = gene_A, to = gene_B),
  directed = FALSE
)

palette_colors <- c(
  "#66C2A5", 
  "#8DA0CB", 
  "#F87262", 
  "#A6D854", 
  "#E78AC3", 
  "#FFA07A",
  "#E5C494", 
  "#B3B3B3",
  "#C6A5E2"
)


p <- ggraph(ppi_graph, layout = "fr") + 
  geom_edge_link(color = "grey60") +
  
  geom_node_point(
    aes(fill = ppi_module, shape = type, color = hub_gene), 
    size = 5,
    stroke = 1
  ) +
  
  geom_node_text(
    aes(label = Gene),  
    size = 4.5,
    repel = TRUE,
    family = "Arial",
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.3, "lines")
  ) +
  
  scale_shape_manual(
    values = c("DEP" = 24, "not-DEP" = 21),
    name = "Protein type"
  ) +
  
  scale_fill_manual(
    values = palette_colors,
    name = "SubModule",
    guide = "none"
  )  +
  
  scale_color_manual(
    values = c("Yes" = "red", "No" = "grey70"),
    name = "Hub Gene"
  ) +
  
  theme_void() +
  
  theme(
    legend.position = 'bottom'
  ) +
  
  labs(
    title = NULL,
    fill = "SubModule"  # 显式指定 fill 图例名
  ) +
  
  guides(
    shape = guide_legend(override.aes = list(fill = "grey70", size = 5))
  )

p

ggsave(paste0(out_dir,'plot/ppi/ppi_all.pdf'),
       plot = p,
       width = 12,
       height = 12,
       device = cairo_pdf)

write.csv(ppi_mod_edge_all,paste0(out_dir, 'table/module/ppi_mod_edge.csv'))
write.csv(ppi_mod_node_all,paste0(out_dir, 'table/module/ppi_mod_node.csv'))


############################# module preservation
# RNA WGCNA
load('/data/projects/China_Brain_MultiOmics/humanBrain_RNAseq/CBMAP_RNAseq_protein_coding/WGCNA/RNA_seq_gene_name_expression.Rdata')
main_expr <- as.data.frame(main_expr)
rna_seq <- main_expr[rownames(main_expr)%in%rownames(pro_raw_50_im),colnames(main_expr)%in%pro_50_raw$gene_name]
colnames(rna_seq) <- pro_50_raw$protein[match(colnames(rna_seq),pro_50_raw$gene_name)]
pro_use <- pro_raw_50_im[rownames(pro_raw_50_im)%in%rownames(rna_seq),colnames(pro_raw_50_im)%in%colnames(rna_seq)]

multiExpr <- list(
  Protein = list(data = pro_use), 
  RNA = list(data = rna_seq)
)

moduleLabels <- adnc_lim_all$net$colors 
moduleColors <- labels2colors(moduleLabels)
names(moduleColors) <- names(moduleLabels)
moduleColors <- moduleColors[colnames(pro_use)]

colorList <- list(Protein = moduleColors)

set.seed(123)
mp_result <- modulePreservation(
  multiData = multiExpr,
  multiColor = colorList,
  referenceNetworks = 1,
  nPermutations = 100,
  networkType = "signed",
  corFnc = "bicor",
  verbose = 3
)
preservationStats <- as.data.frame(mp_result$preservation$Z[[1]][[2]])
preservationStats$module <- rownames(preservationStats)
preservationStats$module_num <- adnc_lim_all$color_table$ModuleNumber[match(preservationStats$module, adnc_lim_all$color_table$ModuleColor)]
write.xlsx(preservationStats, paste0(out_dir,'table/module/module_preservation.xlsx'))
preservationStats <- as.data.frame(fread(paste0(out_dir,'module/module_preservation.txt')))

module_preserve <- preservationStats[!preservationStats$module%in%c('grey','gold'),]
module_preserve$module_num <- paste0('ME', module_preserve$module_num)
module_preserve$module_num <- factor(module_preserve$module_num, levels = paste0('ME',1:(nrow(color_table)-1)))

p_module_pre <- ggplot(module_preserve, aes(x = module_num, y = Zsummary.pres)) +
  geom_point(aes(size = moduleSize, fill = module), 
             color = "black", shape = 21, stroke = 0.8) +  # shape=21支持fill颜色
  scale_fill_identity() +  # 使用实际的颜色值
  scale_size_continuous(range = c(1.96, 10)) +
  geom_hline(yintercept = c(1.96, 10), linetype = "dashed", color = c("darkgreen", "red")) +
  scale_y_continuous(breaks = seq(0, max(module_preserve$Zsummary.pres) + 2, by = 2.5)) +  
  labs(
    x = "Module",
    y = "Z-summary",
    size = "Protein Count"
  ) +
  theme_minimal() +
  theme_bw() +
  theme(
    panel.grid = element_line(color = 'white'),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
  )
p_module_pre
ggsave(paste0(out_dir,'plot/module/module_preservation.pdf'),
       plot = p_module_pre,
       height = 6,
       width = 8,
       dpi = 300)

# rna_module_trait 
targetModules <- module_preserve$module
me_table <- data.frame(sampleid = rownames(rna_seq))

for (mod in targetModules) {
  cat("Processing module:", mod, "\n")
  
  names(adnc_lim_all$moduleColors) <- colnames(pro_raw_50_im)
  module_proteins <- names(adnc_lim_all$moduleColors)[adnc_lim_all$moduleColors == mod]
  rna_select <- rna_seq[,colnames(rna_seq)%in%module_proteins]
  
  MEs <- moduleEigengenes(rna_select, colors = rep(mod, ncol(rna_select)))$eigengenes
  me_table[,mod] <- MEs[match(me_table$sampleid,rownames(MEs)),1]
}

rownames(me_table) <- me_table$sampleid

rna_module_trait <- module_trait_cor(me_table)
rna_module_trait$trait_table$pmarker <- ifelse(rna_module_trait$trait_table$adj.P.Val < 0.05, "*", '')
rna_module_trait$trait_table$module_num <- module_preserve$module_num[match(rna_module_trait$trait_table$me, module_preserve$module)]
rna_module_trait$trait_table$module_num <- factor(rna_module_trait$trait_table$module_num, levels = paste0('ME', 1:(nrow(color_table)-1)))
rna_module_trait$trait_table$trait <- factor(rna_module_trait$trait_table$trait, levels = c('Age','LATE','ARTAG','PART','CVD','LBD','Cerad','Braak','Abeta','ADNC'))

p_rna_trait <- ggplot(rna_module_trait$trait_table, aes(x = module_num, y = trait, fill = -log10(P.Value))) +
  geom_tile(color = "white") +
  geom_text(aes(label = pmarker), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = "#fd8d3c", name = "-log10(P)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, margin = margin(b = 4))
  ) +
  labs(x = NULL, y = NULL, title = NULL)
p_rna_trait
ggsave(paste0(out_dir,'plot/module/module_trait.pdf'),
       plot = p_rna_trait,
       height = 6,
       width = 8,
       dpi = 300)

############################ module connectivity
module_table <- adnc_lim_all$color_table[,c(2,4)]
colnames(module_table) <- c('module_color','module_name')

wgcna_network <- function(data, moduleColors, power = power, net, threshold) {
  # Step 1: 计算全基因邻接矩阵（也可替换成 TOM）
  adjacency_matrix <- adjacency(data, power = power, type = "signed")
  
  # Step 2: 导出 Cytoscape 网络
  cyto_data <- exportNetworkToCytoscape(
    adjacency_matrix,
    threshold = 0.2,
    nodeNames = colnames(data),
    nodeAttr = moduleColors
  )
  
  # Step 3: 整理 edge 和 node
  edges <- cyto_data$edgeData
  nodes <- cyto_data$nodeData
  colnames(edges)[3] <- c("weight")
  colnames(nodes)[1:3] <- c("name", 'no',"module_color")
  
  grey_nodes <- nodes$name[nodes$module_color=='grey']
  edges_no_grey <- edges[-which(edges$fromNode %in% grey_nodes | edges$toNode %in% grey_nodes),]
  
  edges_select <- edges_no_grey[edges_no_grey$weight>threshold,]
  edges_select$gene_a <- pro_50_raw$gene_name[match(edges_select$fromNode, pro_50_raw$protein)]
  edges_select$gene_b <- pro_50_raw$gene_name[match(edges_select$toNode, pro_50_raw$protein)]
  
  select_nodes <- union(edges_select$fromNode, edges_select$toNode)
  nodes_select <- nodes[nodes$name %in% select_nodes,]
  nodes_select$gene <- pro_50_raw$gene_name[match(nodes_select$name, pro_50_raw$protein)]
  nodes_select$hub <- ifelse(nodes_select$gene %in% net$hub_gene$Hub_Gene, 'Hub', 'Not-Hub')
  # nodes_select$module <- module_table$module_color[match(nodes_select$module_color, module_table$module_color)]
  nodes_select$module <- factor(nodes_select$module_color)
  nodes_select$module <- factor(nodes_select$module, levels = adnc_lim_all$color_table$ModuleColor[-1])
  
  graph <- tbl_graph(
    nodes = nodes_select %>%
      dplyr::select(name = gene, module, hub),
    edges = edges_select %>% rename(from = gene_a, to = gene_b, weight = weight),
    directed = FALSE
  )
  
  p <- ggraph(graph, layout = "fr") + 
    geom_edge_link(aes(edge_width = weight), color = "grey50") +
    
    geom_node_point(
      aes(fill = module, shape = hub),
      size = 4, color = "grey", stroke = 0.5
    ) +
    
    geom_node_text(
      aes(label = name),
      size = 3, repel = TRUE, family = "Arial",
      box.padding = unit(0.4, "lines"), point.padding = unit(0.3, "lines")
    ) +
    
    # 使用实际颜色（字符串）作为 fill，不手动指定 color 取值
    scale_fill_identity(name = "Module", guide = "legend") +
    
    scale_shape_manual(
      values = c('Hub' = 24, "Not-Hub" = 21),
      name = "Hub"
    ) +
    
    scale_edge_width_continuous(
      range = c(0.3, 2),
      name = "Weight"
    ) +
    
    guides(
      fill = guide_legend(override.aes = list(shape = 21, color = "black")),
      shape = guide_legend(order = 1),
      edge_width = guide_legend(order = 2)
    ) +
    
    theme_void() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.5, "cm"),
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    )
  
  print(p)
  
  return(
    list(
      edges_select = edges_select,
      nodes_select = nodes_select,
      plot = p
    )
  )
}
# protein
adnc_p_net_wgcn_7 <- wgcna_network(data = pro_raw_50_im,
                                   moduleColors = adnc_lim_all$moduleColors,
                                   power = adnc_lim_all$power,
                                   net = adnc_lim_all,
                                   threshold = 0.7)
p1 <- adnc_p_net_wgcn_7$plot
ggsave(paste0(out_dir,'plot/module/module_connect.png'),
       p1,
       width = 20,
       height = 12,
       dpi = 300)






