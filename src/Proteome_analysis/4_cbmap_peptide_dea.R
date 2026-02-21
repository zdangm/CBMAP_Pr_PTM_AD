library(data.table)
library(readxl)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(writexl)
library(tidyverse)
library(limma)
library(RColorBrewer)
library(impute)
library(biomaRt)
library(openxlsx)
library(readxl)
library(scales)

source('/share/home/sunly/Rscript/cbmap_pro/cbmap_586/2_sample_process.R')
out_dir <- '/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/'

peptide_raw <- fread('/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/pr_jingjie_966_20250214/XB07045B4DA_966samples_Preliminary_analysis_result/peptide_quant.txt')
peptide_raw <- data.frame(peptide_raw[,c(1,3,5,9,999:1985)])
sample_name_pro <- fread("processed_data/human_brain_pro_log_729.txt",nrow=0)
batch_info <- data.frame(fread('processed_data/batch_info.txt'))
sample_adnc$ipb <- batch_info$IP[match(sample_adnc$id, batch_info$sample)]
sample_adnc$enb <- batch_info$enzymolysis[match(sample_adnc$id, batch_info$sample)]

###################### limma peptide adnc
peptide_all <- peptide_raw
colnames(peptide_all)[5:ncol(peptide_all)] <- sub(".*\\.", "", colnames(peptide_all)[5:ncol(peptide_all)])
colnames(peptide_all)[which(colnames(peptide_all)=='repeat2')] <- 'PTB087'
selelct_colnames <- c("Sequence", "Protein.accession", "Gene.name", "Unique..yes.no.", colnames(sample_name_pro)[4:ncol(sample_name_pro)])
peptide_all <- peptide_all[,colnames(peptide_all)%in%selelct_colnames]
peptide_all <- peptide_all[!duplicated(peptide_all$Sequence),]
# peptide_intensity <- peptide_all
peptide_all[, 5:ncol(peptide_all)] <- log2(peptide_all[, 5:ncol(peptide_all)] + 1)
peptide_all$missing_rate <- apply(peptide_all[,5:ncol(peptide_all)], 1, function(x) mean(is.na(x)))
peptide_all <- peptide_all[peptide_all$missing_rate<0.5,]
rownames(peptide_all) <- peptide_all$Sequence

# app_intensity <- data.frame(t(peptide_intensity[peptide_intensity$Sequence=='LVFFAEDVGSNK',c(5:ncol(peptide_intensity))]))
# app_log <- data.frame(t(peptide_all[peptide_all$Sequence=='LVFFAEDVGSNK',c(5:733)]))
# app_rint <- data.frame(t(peptide_all_im[peptide_all_im$Sequence=='LVFFAEDVGSNK',c(5:733)]))
# 
# peptide_app <- data.frame(sample = colnames(peptide_intensity)[5:733])
# peptide_app$intensity <- app_intensity$X33976[match(peptide_app$sample, rownames(app_intensity))]
# peptide_app$log <- app_log$LVFFAEDVGSNK[match(peptide_app$sample, rownames(app_log))]
# peptide_app$rint <- app_rint$LVFFAEDVGSNK[match(peptide_app$sample, rownames(app_rint))]
# write.table(peptide_app, 'processed_data/tmp/APP_LVFFAEDVGSNK.txt', quote = F, row.names = F, sep = '\t')

set.seed(20250514)
peptide_all_im <- impute.knn(as.matrix(peptide_all[,5:(ncol(peptide_all)-1)]), k = 10, rng.seed = 2025)$data %>% as.data.frame()
peptide_all_im <- sweep(peptide_all_im, 2, apply(peptide_all_im, 2, median, na.rm = TRUE), FUN = "-")
rint_transform <- function(x) {
  n <- length(x)
  ranks <- rank(x, ties.method = "average")
  qnorm((ranks - 0.5) / n)
}
peptide_all_im <- apply(peptide_all_im, 1, rint_transform)
peptide_all_im <- t(peptide_all_im)
peptide_all_im <- cbind(peptide = rownames(peptide_all_im), as.data.frame(peptide_all_im))

peptide_all_im <- merge(peptide_all[,1:4], peptide_all_im, by.x = 'Sequence', by.y = 'peptide', all.y = TRUE, all.x = FALSE)
rownames(peptide_all_im) <- peptide_all_im$Sequence
write.table(peptide_all_im, 'processed_data/human_brain_peptide_rint.txt', sep = '\t', quote = F,row.names = F)
peptide_all_im <- as.data.frame(fread('processed_data/human_brain_peptide_rint.txt'))
rownames(peptide_all_im) <- peptide_all_im$Sequence

adnc_id <- intersect(colnames(peptide_all_im), sample_adnc$id)
peptide_all_adnc <- peptide_all_im[,colnames(peptide_all_im)%in%adnc_id]
peptide_all_adnc <- peptide_all_adnc[,match(adnc_id, colnames(peptide_all_adnc))]
sample_adnc <- sample_adnc[sample_adnc$id %in% adnc_id,]
sample_adnc <- sample_adnc[match(adnc_id, sample_adnc$id),]

design <- model.matrix( ~ adnc_num + age + sex_male + PMD + RIN, data = sample_adnc)
fit <- lmFit(peptide_all_adnc, design)
fit <- eBayes(fit)
peptide_all_adnc_limma_re <- topTable(fit, coef = "adnc_num", number = Inf, adjust.method = "BH")
peptide_all_adnc_limma_re$peptide <- rownames(peptide_all_adnc_limma_re)
peptide_all_adnc_limma_re$gene_name <- peptide_all$Gene.name[match(peptide_all_adnc_limma_re$peptide, peptide_all$Sequence)]
peptide_all_adnc_limma_re$pro <- peptide_all$Protein.accession[match(peptide_all_adnc_limma_re$peptide, peptide_all$Sequence)]

length(peptide_all_adnc_limma_re$peptide[peptide_all_adnc_limma_re$adj.P.Val<0.05&peptide_all_adnc_limma_re$logFC<0])
length(unique(peptide_all_adnc_limma_re$gene_name[peptide_all_adnc_limma_re$adj.P.Val<0.05]))

write.xlsx(peptide_all_adnc_limma_re, paste0(out_dir, 'table/dea_result/cbmap_adnc_lim_peptide.xlsx'))


########################## limma peptide adjust cov 
# adnc
sample_adnc_cov <- sample_adnc
sample_adnc_cov$lbd[is.na(sample_adnc_cov$lbd)] <- 0
sample_adnc_cov$cvd[is.na(sample_adnc_cov$cvd)] <- 0
sample_adnc_cov$late[is.na(sample_adnc_cov$late)] <- 0
adnc_id_cov <- intersect(sample_adnc_cov$id, colnames(pro_50_raw)[4:ncol(pro_50_raw)])
peptide_adnc_limma_cov <- peptide_all_im[,colnames(peptide_all_im)%in%adnc_id_cov]
peptide_adnc_limma_cov <- peptide_adnc_limma_cov[, match(adnc_id_cov, colnames(peptide_adnc_limma_cov))]
sample_adnc_cov <- sample_adnc_cov[sample_adnc_cov$id%in%adnc_id_cov,]
sample_adnc_cov <- sample_adnc_cov[match(adnc_id_cov, sample_adnc_cov$id), ]
summary(sample_adnc_cov$adnc_fac)

design <- model.matrix( ~ adnc_num + age + sex_male + PMD + RIN + lbd + cvd + late, data = sample_adnc_cov)
fit <- lmFit(peptide_adnc_limma_cov, design)
fit <- eBayes(fit)
adnc_num_lim_re_cov <- topTable(fit, coef = "adnc_num", number = Inf, adjust.method = "BH")
adnc_num_lim_re_cov$peptide <- rownames(adnc_num_lim_re_cov)
adnc_num_lim_re_cov$gene_name <- peptide_all_im$Gene.name[match(adnc_num_lim_re_cov$peptide, peptide_all_im$Sequence)] # PLCB1过fdr校正，P值均小于0.05
adnc_num_lim_re_cov$pro <- peptide_all_im$Protein.accession[match(adnc_num_lim_re_cov$peptide, peptide_all_im$Sequence)] # PLCB1过fdr校正，P值均小于0.05

nrow(adnc_num_lim_re_cov[adnc_num_lim_re_cov$adj.P.Val<0.05&adnc_num_lim_re_cov$logFC<0,])
length(unique(adnc_num_lim_re_cov$gene_name[adnc_num_lim_re_cov$adj.P.Val<0.05&adnc_num_lim_re_cov$logFC<0]))

write.xlsx(adnc_num_lim_re_cov, paste0(out_dir, 'table/dea_result/cbmap_adnc_lim_peptide_cov.xlsx'))


# braak
braak_id <- intersect(colnames(peptide_all_im), sample_braak$id)
peptide_all_braak <- peptide_all_im[,colnames(peptide_all_im)%in%braak_id]
peptide_all_braak <- peptide_all_braak[,match(braak_id, colnames(peptide_all_braak))]
sample_braak <- sample_braak[sample_braak$id %in% braak_id,]
sample_braak <- sample_braak[match(braak_id, sample_braak$id),]

design <- model.matrix( ~ Braak.NFT.stage + age + sex_male + PMD + RIN, data = sample_braak)
fit <- lmFit(peptide_all_braak, design)
fit <- eBayes(fit)
peptide_all_braak_limma_re <- topTable(fit, coef = "Braak.NFT.stage", number = Inf, adjust.method = "BH")
peptide_all_braak_limma_re$peptide <- rownames(peptide_all_braak_limma_re)
peptide_all_braak_limma_re$gene_name <- peptide_all$Gene.name[match(peptide_all_braak_limma_re$peptide, peptide_all$Sequence)]

length(peptide_all_braak_limma_re$peptide[peptide_all_braak_limma_re$adj.P.Val<0.05&peptide_all_braak_limma_re$logFC>0])
length(unique(peptide_all_braak_limma_re$gene_name[peptide_all_braak_limma_re$adj.P.Val<0.05&peptide_all_braak_limma_re$logFC<0]))

write.xlsx(peptide_all_braak_limma_re, paste0(out_dir, 'table/dea_result/cbmap_braak_lim_peptide.xlsx'))


#abeta
abeta_id <- intersect(colnames(peptide_all_im), sample_abeta$id)
peptide_all_abeta <- peptide_all_im[,colnames(peptide_all_im)%in%abeta_id]
peptide_all_abeta <- peptide_all_abeta[,match(abeta_id, colnames(peptide_all_abeta))]
sample_abeta <- sample_abeta[sample_abeta$id %in% abeta_id,]
sample_abeta <- sample_abeta[match(abeta_id, sample_abeta$id),]

design <- model.matrix( ~ A_beta_0_3 + age + sex_male + PMD + RIN, data = sample_abeta)
fit <- lmFit(peptide_all_abeta, design)
fit <- eBayes(fit)
peptide_all_abeta_limma_re <- topTable(fit, coef = "A_beta_0_3", number = Inf, adjust.method = "BH")
peptide_all_abeta_limma_re$peptide <- rownames(peptide_all_abeta_limma_re)
peptide_all_abeta_limma_re$gene_name <- peptide_all$Gene.name[match(peptide_all_abeta_limma_re$peptide, peptide_all$Sequence)]

length(peptide_all_abeta_limma_re$peptide[peptide_all_abeta_limma_re$adj.P.Val<0.05])
length(peptide_all_abeta_limma_re$peptide[peptide_all_abeta_limma_re$adj.P.Val<0.05&peptide_all_abeta_limma_re$logFC<0])
length(unique(peptide_all_abeta_limma_re$gene_name[peptide_all_abeta_limma_re$adj.P.Val<0.05&peptide_all_abeta_limma_re$logFC>0]))
length(unique(peptide_all_abeta_limma_re$gene_name[peptide_all_abeta_limma_re$adj.P.Val<0.05]))

write.xlsx(peptide_all_abeta_limma_re, paste0(out_dir, 'table/dea_result/cbmap_abeta_lim_peptide.xlsx'))


########################## sensitivity analysis
# adnc
sample_adnc_sen <- sample_adnc[sample_adnc$bank%in%c('zju','pumc'),]
adnc_id_sen <- intersect(colnames(peptide_all_im), sample_adnc_sen$id)
peptide_all_sen_adnc <- peptide_all_im[,colnames(peptide_all_im)%in%adnc_id_sen]
peptide_all_sen_adnc <- peptide_all_sen_adnc[,match(adnc_id_sen, colnames(peptide_all_sen_adnc))]
sample_adnc_sen <- sample_adnc_sen[sample_adnc_sen$id %in% adnc_id_sen,]
sample_adnc_sen <- sample_adnc_sen[match(adnc_id_sen, sample_adnc_sen$id),]

design <- model.matrix( ~ adnc_num + age + sex_male + PMD + RIN, data = sample_adnc_sen)
fit <- lmFit(peptide_all_sen_adnc, design)
fit <- eBayes(fit)
peptide_all_adnc_limma_sen_re <- topTable(fit, coef = "adnc_num", number = Inf, adjust.method = "BH")
peptide_all_adnc_limma_sen_re$peptide <- rownames(peptide_all_adnc_limma_sen_re)
peptide_all_adnc_limma_sen_re$gene_name <- peptide_all$Gene.name[match(peptide_all_adnc_limma_sen_re$peptide, peptide_all$Sequence)]
peptide_all_adnc_limma_sen_re$pro <- peptide_all$Protein.accession[match(peptide_all_adnc_limma_sen_re$peptide, peptide_all$Sequence)]

length(peptide_all_adnc_limma_sen_re$peptide[peptide_all_adnc_limma_sen_re$adj.P.Val<0.05&peptide_all_adnc_limma_sen_re$logFC<0])
length(unique(peptide_all_adnc_limma_sen_re$gene_name[peptide_all_adnc_limma_sen_re$adj.P.Val<0.05]))

write.xlsx(peptide_all_adnc_limma_sen_re, paste0(out_dir, 'table/dea_result/cbmap_adnc_lim_peptide_sen.xlsx'))

# braak
sample_braak_sen <- sample_braak[sample_braak$bank%in%c('zju','pumc'),]
braak_id_sen <- intersect(colnames(peptide_all_im), sample_braak_sen$id)
peptide_all_sen_braak <- peptide_all_im[,colnames(peptide_all_im)%in%braak_id_sen]
peptide_all_sen_braak <- peptide_all_sen_braak[,match(braak_id_sen, colnames(peptide_all_sen_braak))]
sample_braak_sen <- sample_braak_sen[sample_braak_sen$id %in% braak_id_sen,]
sample_braak_sen <- sample_braak_sen[match(braak_id_sen, sample_braak_sen$id),]

design <- model.matrix( ~ Braak.NFT.stage + age + sex_male + PMD + RIN, data = sample_braak_sen)
fit <- lmFit(peptide_all_sen_braak, design)
fit <- eBayes(fit)
peptide_all_braak_limma_sen_re <- topTable(fit, coef = "Braak.NFT.stage", number = Inf, adjust.method = "BH")
peptide_all_braak_limma_sen_re$peptide <- rownames(peptide_all_braak_limma_sen_re)
peptide_all_braak_limma_sen_re$gene_name <- peptide_all$Gene.name[match(peptide_all_braak_limma_sen_re$peptide, peptide_all$Sequence)]
peptide_all_braak_limma_sen_re$pro <- peptide_all$Protein.accession[match(peptide_all_braak_limma_sen_re$peptide, peptide_all$Sequence)]

length(peptide_all_braak_limma_sen_re$peptide[peptide_all_braak_limma_sen_re$adj.P.Val<0.05&peptide_all_braak_limma_sen_re$logFC<0])
length(unique(peptide_all_braak_limma_sen_re$gene_name[peptide_all_braak_limma_sen_re$adj.P.Val<0.05]))

write.xlsx(peptide_all_braak_limma_sen_re, paste0(out_dir, 'table/dea_result/cbmap_braak_lim_peptide_sen.xlsx'))

# abeta
sample_abeta_sen <- sample_abeta[sample_abeta$bank%in%c('zju','pumc'),]
abeta_id_sen <- intersect(colnames(peptide_all_im), sample_abeta_sen$id)
peptide_all_sen_abeta <- peptide_all_im[,colnames(peptide_all_im)%in%abeta_id_sen]
peptide_all_sen_abeta <- peptide_all_sen_abeta[,match(abeta_id_sen, colnames(peptide_all_sen_abeta))]
sample_abeta_sen <- sample_abeta_sen[sample_abeta_sen$id %in% abeta_id_sen,]
sample_abeta_sen <- sample_abeta_sen[match(abeta_id_sen, sample_abeta_sen$id),]

design <- model.matrix( ~ A_beta_0_3 + age + sex_male + PMD + RIN, data = sample_abeta_sen)
fit <- lmFit(peptide_all_sen_abeta, design)
fit <- eBayes(fit)
peptide_all_abeta_limma_sen_re <- topTable(fit, coef = "A_beta_0_3", number = Inf, adjust.method = "BH")
peptide_all_abeta_limma_sen_re$peptide <- rownames(peptide_all_abeta_limma_sen_re)
peptide_all_abeta_limma_sen_re$gene_name <- peptide_all$Gene.name[match(peptide_all_abeta_limma_sen_re$peptide, peptide_all$Sequence)]
peptide_all_abeta_limma_sen_re$pro <- peptide_all$Protein.accession[match(peptide_all_abeta_limma_sen_re$peptide, peptide_all$Sequence)]

length(peptide_all_abeta_limma_sen_re$peptide[peptide_all_abeta_limma_sen_re$adj.P.Val<0.05&peptide_all_abeta_limma_sen_re$logFC<0])
length(unique(peptide_all_abeta_limma_sen_re$gene_name[peptide_all_abeta_limma_sen_re$adj.P.Val<0.05]))

write.xlsx(peptide_all_abeta_limma_sen_re, paste0(out_dir, 'table/dea_result/cbmap_abeta_lim_peptide_sen.xlsx'))


########################## miami plot
pathway_func <- read_xlsx('splicing_pathway_categories.xlsx')
sheet_names <- excel_sheets(paste0(out_dir, 'table/enrich_result/cbmap_adnc_gkg.xlsx'))
enrichment <- lapply(sheet_names, function(sheet) {
  read_excel(paste0(out_dir, 'table/enrich_result/cbmap_adnc_gkg.xlsx'),sheet = sheet)
})

adnc_go <- enrichment[[3]]
adnc_go$Category <- pathway_func$Category[match(adnc_go$Description, pathway_func$Term)]

category <- unique(adnc_go$Category)
gene_categray <- adnc_go %>%
  separate_rows(Gene_Symbols, sep = ', ') %>%
  distinct(Category, Gene_Symbols)  # 去重
gene_categray$Category <- ifelse(gene_categray$Category=="Cytoskeleton & contractile structure", "Cytoskeleton", gene_categray$Category)
gene_categray$cat <- ifelse(gene_categray$Category%in%c("RNA processing & splicing","Synapse & neuron","Cytoskeleton","Cell adhesion & signaling"), gene_categray$Category, NA)


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes <- unique(peptide_all_adnc_limma_re$gene_name)
gene_anno <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = genes,
  mart = ensembl
)

valid_chr <- c(as.character(1:22), "X", "Y")
gene_anno <- gene_anno[gene_anno$chromosome_name %in% valid_chr, ]
gene_anno$chromosome_name <- factor(gene_anno$chromosome_name, levels = valid_chr)

dual_miami_plot <- function(data, gene_anno, fdr_col = "P.Value", logfc_col = "logFC", genename_col = "gene_name",
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
  data$label <- paste0(data$gene_name, '_', data$peptide)
  
  # 加入富集功能
  gene_cat <- gene_categray[!is.na(gene_categray$cat),]
  
  priority <- c("RNA processing & splicing", "Cytoskeleton", "Synapse & neuron", "Cell adhesion & signaling")
  gene_cat <- gene_cat %>%
    dplyr::mutate(Category = factor(Category, levels = priority)) %>%
    dplyr::arrange(Gene_Symbols, Category) %>%
    dplyr::group_by(Gene_Symbols) %>%
    dplyr::slice(1)
  
  gene_cat$genename <- gene_cat$Gene_Symbols
  
  data <- data %>%
    left_join(gene_cat[,c(1,4)], by = c("gene_name" = "genename"))
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
    "RNA processing & splicing" = "#C3D9F0",         # 淡蓝
    "Cytoskeleton" = "#FFF2C2", # 淡黄
    "Synapse & neuron" = "#B5E6DE",                 # 淡青绿
    "Cell adhesion & signaling" = "#D3C9EB"         # 淡紫
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
      aes(x = pos_cum, y = logp_dir, label = label, fill = Category),
      color = 'black',
      size = 4.5,
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
      aes(x = pos_cum, y = logp_dir, label = label, fill = 'white'),
      color = 'black',
      size = 4.5,
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
    # scale_color_manual(values = c(category_colors, chr_colors_up), na.translate = FALSE, guide = "none") +
    
    scale_y_continuous(name = "Positive -log10(P)", limits = c(0, maxp)) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    coord_cartesian(xlim = x_limits) +
    
    theme_classic() +
    theme(
      # legend.position = "right",
      # plot.margin = margin(t = 5, r = 5, b = 5, l = 10),
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
      aes(x = pos_cum, y = logp_dir, label = label, fill = Category),
      color = 'black',
      size = 4.5,
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
      aes(x = pos_cum, y = logp_dir, label = label, fill = 'white'),
      color = 'black',
      size = 4.5,
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
    scale_fill_manual(values = category_colors, 
                      na.translate = FALSE, 
                      name = NULL, 
                      drop = FALSE,
                      guide = guide_legend( # 同样修改
                        override.aes = list(label = "", size = 0),
                        keywidth = unit(12, "pt"),
                        keyheight = unit(12, "pt")
                      )) +
    # scale_color_manual(values = c(category_colors, chr_colors_down), na.translate = FALSE, guide = "none") +
    
    scale_y_continuous(name = "Negative -log10(P)", limits = c(-maxp,0),labels = abs) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, position = 'top') +
    coord_cartesian(xlim = x_limits) +
    
    theme_classic() +
    theme(
      # legend.position = "bottom",
      # legend.box = "horizontal",
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
  combined_plot <- combined_plot & theme(legend.position = "bottom", legend.box = "horizontal")
  print(combined_plot)
  return(combined_plot)
}

adnc_plot <- dual_miami_plot(data = peptide_all_adnc_limma_re, gene_anno = gene_anno)
ggsave(paste0(out_dir, 'plot/miami/adnc_peptide.pdf'),
       adnc_plot,
       width = 12,
       height = 9)
 
adnc_cov_plot <- dual_miami_plot(data = adnc_num_lim_re_cov, gene_anno = gene_anno)
ggsave(paste0(out_dir,'plot/miami/adnc_peptide_cov.png'),
       adnc_cov_plot,
       width = 12,
       height = 8,
       dpi = 300)

adnc_sen_plot <- dual_miami_plot(data = peptide_all_adnc_limma_sen_re, gene_anno = gene_anno)
ggsave(paste0(out_dir,'plot/miami/adnc_peptide_sen.png'),
       adnc_sen_plot,
       width = 12,
       height = 8,
       dpi = 300)

braak_plot <- dual_miami_plot(data = peptide_all_braak_limma_re, gene_anno = gene_anno)
ggsave(paste0(out_dir, 'plot/miami/braak_peptide.png'),
       braak_plot,
       width = 12,
       height = 8,
       dpi = 300)

braak_sen_plot <- dual_miami_plot(data = peptide_all_braak_limma_sen_re, gene_anno = gene_anno)
ggsave(paste0(out_dir,'plot/miami/braak_peptide_sen.pdf'),
       braak_sen_plot,
       width = 12,
       height = 8)

abeta_plot <- dual_miami_plot(data = peptide_all_abeta_limma_re, gene_anno = gene_anno)
ggsave(paste0(out_dir,'plot/miami/abeta_peptide.png'),
       abeta_plot,
       width = 12,
       height = 8,
       dpi = 300)

abeta_sen_plot <- dual_miami_plot(data = peptide_all_abeta_limma_sen_re, gene_anno = gene_anno)
ggsave(paste0(out_dir,'plot/miami/abeta_peptide_sen.pdf'),
       abeta_sen_plot,
       width = 12,
       height = 8)

################### top proteins box plot
peptide_all_adnc_limma_re$gene_peptide <- paste0(peptide_all_adnc_limma_re$gene_name,'_',peptide_all_adnc_limma_re$peptide)
select_peptide <- c('LVFFAEDVGSNK','FSEPDPSHTLEER','EVSSATNALR','EFEVYGPIK')

selece_peptide_intensity <- peptide_all[rownames(peptide_all)%in%select_peptide,colnames(peptide_all)%in%sample_adnc$id]
selece_peptide_intensity <- data.frame(t(selece_peptide_intensity))
selece_peptide_intensity$id <- rownames(selece_peptide_intensity)
selece_peptide_intensity$adnc_num <- sample_adnc$adnc_num[match(selece_peptide_intensity$id, sample_adnc$id)]
selece_peptide_intensity$bank <- sample_adnc$bank[match(selece_peptide_intensity$id, sample_adnc$id)]

# adnc
for (pro in select_peptide) {
  data <- selece_peptide_intensity[,c('id','adnc_num','bank',pro)]
  data$Disease <- ifelse(data$adnc_num==0,'HC',
                         ifelse(data$adnc_num==1,'ADNC_L',
                                ifelse(data$adnc_num==2,'ADNC_M','ADNC_H')))
  data$Disease <- factor(data$Disease, levels = c('HC','ADNC_L','ADNC_M','ADNC_H'))
  gene_peptide <- peptide_all_adnc_limma_re$gene_peptide[peptide_all_adnc_limma_re$peptide==pro]
  
  y_min <- min(data[, pro], na.rm = TRUE)
  y_max <- max(data[, pro], na.rm = TRUE)
  
  my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(unique(data$Disease))) 
  
  p <- ggplot(data, aes(x = Disease, y = .data[[names(data)[4]]], fill = Disease)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene_peptide, '_adnc'), x = "", y = "")
  print(p)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_adnc_all.pdf'), plot = p, width = 5, height = 5)
  
  p1 <- ggplot(data[data$bank %in% c('zju','csu'),], aes(x = Disease, y = .data[[names(data)[4]]], fill = Disease)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title =paste0(gene_peptide,'_adnc_zju_csu'), x = "", y = "")
  print(p1)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_adnc_zc.pdf'), plot = p1, width = 5, height = 5)
  
  p2 <- ggplot(data[data$bank=='pumc',], aes(x = Disease, y = .data[[names(data)[4]]], fill = Disease)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene_peptide,'_adnc_pumc'), x = "", y = "")
  print(p2)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_adnc_pumc.pdf'), plot = p2, width = 5, height = 5)
}

# braak
peptide_all_braak_limma_re$gene_peptide <- paste0(peptide_all_braak_limma_re$gene_name,'_',peptide_all_braak_limma_re$peptide)
select_peptide <- c('LVFFAEDVGSNK','FTDYCDLNK','EVSSATNALR','VNYDTTESK')

selece_peptide_intensity <- peptide_all[rownames(peptide_all)%in%select_peptide,colnames(peptide_all)%in%sample_braak$id]
selece_peptide_intensity <- data.frame(t(selece_peptide_intensity))
selece_peptide_intensity$id <- rownames(selece_peptide_intensity)
selece_peptide_intensity$braak <- sample_braak$Braak.NFT.stage[match(selece_peptide_intensity$id, sample_braak$id)]
selece_peptide_intensity$bank <- sample_braak$bank[match(selece_peptide_intensity$id, sample_braak$id)]
selece_peptide_intensity$braak <- factor(selece_peptide_intensity$braak)
selece_peptide_intensity$braak_type <- ifelse(selece_peptide_intensity$braak==0,'B0',
                                              ifelse(selece_peptide_intensity$braak%in%c(1,2),'B1',
                                                     ifelse(selece_peptide_intensity$braak%in%c(3,4),'B2','B3')))
selece_peptide_intensity$braak_type <- factor(selece_peptide_intensity$braak_type, levels = c('B0','B1','B2','B3'))

for (i in select_peptide) {
  data <- selece_peptide_intensity[,c('id','braak','braak_type','bank',i)]
  gene_peptide <- peptide_all_braak_limma_re$gene_peptide[peptide_all_braak_limma_re$peptide==i]
  
  y_min <- min(data[, i], na.rm = TRUE)
  y_max <- max(data[, i], na.rm = TRUE)
  
  my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(unique(data$braak_type)))
  
  p <- ggplot(data, aes(x = braak_type, y = .data[[names(data)[5]]], fill = braak_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene_peptide, '_braak'), x = "", y = "")
  print(p)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_braak_all.pdf'), plot = p, width = 5, height = 5)
  
  p1 <- ggplot(data[data$bank %in% c('zju','csu'),], aes(x = braak_type, y = .data[[names(data)[5]]], fill = braak_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title =paste0(gene_peptide,'_braak_zju_csu'), x = "", y = "")
  print(p1)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_braak_zc.pdf'), plot = p1, width = 5, height = 5)
  
  p2 <- ggplot(data[data$bank=='pumc',], aes(x = braak_type, y = .data[[names(data)[5]]], fill = braak_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene_peptide,'_braak_pumc'), x = "", y = "")
  print(p2)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_braak_pumc.pdf'), plot = p2, width = 5, height = 5)
  
}

# abeta
peptide_all_abeta_limma_re$gene_peptide <- paste0(peptide_all_abeta_limma_re$gene_name,'_',peptide_all_abeta_limma_re$peptide)
select_peptide <- c('LVFFAEDVGSNK','FTDYCDLNK','EVSSATNALR','EFEVYGPIK')

selece_peptide_intensity <- peptide_all[rownames(peptide_all)%in%select_peptide,colnames(peptide_all)%in%sample_abeta$id]
selece_peptide_intensity <- data.frame(t(selece_peptide_intensity))
selece_peptide_intensity$id <- rownames(selece_peptide_intensity)
selece_peptide_intensity$abeta <- sample_abeta$A_beta_0_3[match(selece_peptide_intensity$id, sample_abeta$id)]
selece_peptide_intensity$bank <- sample_abeta$bank[match(selece_peptide_intensity$id, sample_abeta$id)]
selece_peptide_intensity$abeta <- factor(selece_peptide_intensity$abeta)
selece_peptide_intensity$abeta_type <- ifelse(selece_peptide_intensity$abeta==0,'A0',
                                              ifelse(selece_peptide_intensity$abeta==1,'A1',
                                                     ifelse(selece_peptide_intensity$abeta==2,'A2','A3')))
selece_peptide_intensity$abeta_type <- factor(selece_peptide_intensity$abeta_type, levels = c('A0','A1','A2','A3'))

for (i in select_peptide) {
  data <- selece_peptide_intensity[,c('id','abeta','abeta_type','bank',i)]
  gene_peptide <- peptide_all_abeta_limma_re$gene_peptide[peptide_all_abeta_limma_re$peptide==i]
  
  y_min <- min(data[, i], na.rm = TRUE)
  y_max <- max(data[, i], na.rm = TRUE)
  
  my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(unique(data$abeta)))
  
  p <- ggplot(data, aes(x = abeta_type, y = .data[[names(data)[5]]], fill = abeta_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene_peptide, '_abeta'), x = "", y = "")
  print(p)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_abeta_all.pdf'), plot = p, width = 5, height = 5)
  
  p1 <- ggplot(data[data$bank %in% c('zju','csu'),], aes(x = abeta_type, y = .data[[names(data)[5]]], fill = abeta_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title =paste0(gene_peptide,'_abeta_zju_csu'), x = "", y = "")
  print(p1)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_abeta_zc.pdf'), plot = p1, width = 5, height = 5)
  
  p2 <- ggplot(data[data$bank=='pumc',], aes(x = abeta_type, y = .data[[names(data)[5]]], fill = abeta_type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
    geom_jitter(width = 0.2, alpha = 0.6, size = 1.1, color = "black") +  
    scale_fill_manual(values = my_colors) +  
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    coord_cartesian(ylim = c(y_min, y_max)) + 
    theme_minimal() +
    theme_bw() +
    theme(panel.grid = element_line(color = 'white'),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          plot.title = element_text(size = 16),
          legend.position = "none") +  
    labs(title = paste0(gene_peptide,'_abeta_pumc'), x = "", y = "")
  print(p2)
  ggsave(paste0(out_dir,'plot/box_plot/',gene_peptide,'_abeta_pumc.pdf'), plot = p2, width = 5, height = 5)
}


###app_LVFFAEDVGSNK双峰分布的原因
# adnc bank
peptide_APP_limma_t <- data.frame(t(peptide_APP_limma))
peptide_APP_limma_t$id <- rownames(peptide_APP_limma_t)

# bank
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$bank=='zju']])
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$bank=='csu']])
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$bank=='pumc']]) # pumc有明显的双峰分布
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$bank%in%c('zju','csu')]])
# age

# adnc
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$adnc_num==0]],
     main = 'APP_LVFFAEDVGSNK_ADNC_HC',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$adnc_num==1]],
     main = 'APP_LVFFAEDVGSNK_ADNC_L',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$adnc_num==2]],
     main = 'APP_LVFFAEDVGSNK_ADNC_M',
     xlab = '',
     ylab = '') 
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$adnc_num==3]],
     main = 'APP_LVFFAEDVGSNK_ADNC_H',
     xlab = '',
     ylab = '')

hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$adnc_num%in%c(0,1)]]) 
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$adnc_num%in%c(2,3)]])
hist(peptide_APP_limma_t$LVFFAEDVGSNK[rownames(peptide_APP_limma_t)%in%sample_adnc$id[sample_adnc$adnc_num%in%c(0,1,2)]])

# age、pmd、rin
library(mclust)
library(impute)
# peptide_APP_limma_t_im <- impute.knn(as.matrix(peptide_APP_limma_t), k = 10, rng.seed= 2025)$data %>% as.data.frame()
peptide_APP_limma_t_im <- peptide_APP_limma_t[,c('LVFFAEDVGSNK','id')]
peptide_APP_limma_t_im <- na.omit(peptide_APP_limma_t_im)
# model <- Mclust(peptide_APP_limma_t_im$LVFFAEDVGSNK)
# plot(model, what = "density")
model <- Mclust(peptide_APP_limma_t_im$LVFFAEDVGSNK, G = 2)
peptide_APP_limma_t_im$cluster <- as.factor(model$classification)
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK,
     main = 'APP_LVFFAEDVGSNK_all',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$cluster==1],
     main = 'APP_LVFFAEDVGSNK_cluster1',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$cluster==2],
     main = 'APP_LVFFAEDVGSNK_cluster2',
     xlab = '',
     ylab = '')


# peptide_APP_limma_t_im$id <- rownames(peptide_APP_limma_t_im)
peptide_APP_limma_t_im$age <- sample_adnc$age[match(peptide_APP_limma_t_im$id, sample_adnc$id)]
peptide_APP_limma_t_im$age_group <- ifelse(peptide_APP_limma_t_im$age<60,'<60',
                                           ifelse(peptide_APP_limma_t_im$age>=60&peptide_APP_limma_t_im$age<70,'60-70',
                                                  ifelse(peptide_APP_limma_t_im$age>=70&peptide_APP_limma_t_im$age<80,'70-80',
                                                         ifelse(peptide_APP_limma_t_im$age>=80&peptide_APP_limma_t_im$age<90,'80-90','>90'))))
peptide_APP_limma_t_im$age_group <- factor(peptide_APP_limma_t_im$age_group, levels = c('<60','60-70','70-80','80-90','>90'))
peptide_APP_limma_t_im$bank <- sample_adnc$bank[match(peptide_APP_limma_t_im$id, sample_adnc$id)]
peptide_APP_limma_t_im$bank <- factor(peptide_APP_limma_t_im$bank)
peptide_APP_limma_t_im$pmd <- sample_adnc$PMD[match(peptide_APP_limma_t_im$id, sample_adnc$id)]
peptide_APP_limma_t_im$rin <- sample_adnc$RIN[match(peptide_APP_limma_t_im$id, sample_adnc$id)]
peptide_APP_limma_t_im$ipb <- sample_adnc$ipb[match(peptide_APP_limma_t_im$id, sample_adnc$id)]
peptide_APP_limma_t_im$ipb <- factor(peptide_APP_limma_t_im$ipb, levels = c(paste0('B',1:21)))
peptide_APP_limma_t_im$enb <- sample_adnc$enb[match(peptide_APP_limma_t_im$id, sample_adnc$id)]
peptide_APP_limma_t_im$enb <- factor(peptide_APP_limma_t_im$enb, levels = c(paste0('A',1:4)))
peptide_APP_limma_t_im$adnc <- sample_adnc$adnc_num[match(peptide_APP_limma_t_im$id, sample_adnc$id)]
peptide_APP_limma_t_im$adnc_fac <- factor(peptide_APP_limma_t_im$adnc, levels = c(0,1,2,3))
peptide_APP_limma_t_im$sex <- sample_adnc$sex_male[match(peptide_APP_limma_t_im$id, sample_adnc$id)]


ggplot(peptide_APP_limma_t_im, aes(x = age, fill = cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(peptide_APP_limma_t_im, aes(x = sex, fill = cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(peptide_APP_limma_t_im, aes(x = rin, fill = cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(peptide_APP_limma_t_im, aes(x = pmd, fill = cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(peptide_APP_limma_t_im, aes(x = bank, fill = cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

ggplot(peptide_APP_limma_t_im, aes(x = adnc, fill = cluster)) +
  geom_density(alpha = 0.5) +
  theme_minimal()

# 多因素logistic回归，adnc与双峰分布有关
glm_result <- glm(cluster ~ age + adnc_fac + sex + rin + pmd + bank + ipb, data = peptide_APP_limma_t_im, family = binomial)
summary_result <- summary(glm_result)
coef_summary_fac <- as.data.frame(summary_result$coefficients)
coef_summary_fac$fdr <- p.adjust(coef_summary_fac$`Pr(>|z|)`, method = 'BH')

glm_result <- glm(cluster ~ age + adnc + sex + rin + pmd + bank + ipb, data = peptide_APP_limma_t_im, family = binomial)
summary_result <- summary(glm_result)
coef_summary <- as.data.frame(summary_result$coefficients)
coef_summary$fdr <- p.adjust(coef_summary$`Pr(>|z|)`, method = 'BH')

glm_result <- glm(cluster ~ age + sex + rin + pmd + bank + ipb, data = peptide_APP_limma_t_im, family = binomial)
summary_result <- summary(glm_result)
coef_summary_age <- as.data.frame(summary_result$coefficients)
coef_summary_age$fdr <- p.adjust(coef_summary_age$`Pr(>|z|)`, method = 'BH')

ggplot(peptide_APP_limma_t_im, aes(x = adnc, fill = as.factor(cluster))) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", fill = "Cluster") +
  theme_minimal()

ggplot(peptide_APP_limma_t_im, aes(x = age_group, fill = as.factor(cluster))) +
  geom_bar(position = "fill") +
  labs(y = "Proportion", fill = "Cluster") +
  theme_minimal()

####伪时间轨迹
# 控制疾病状态
hc_app <- peptide_APP_limma_t_im[peptide_APP_limma_t_im$adnc==0,]
adnc_app <- peptide_APP_limma_t_im[peptide_APP_limma_t_im$adnc%in%c(1,2,3),]

ggplot(hc_app, aes(x = age, y = LVFFAEDVGSNK)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  # ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_HC", x = "age", y = "intensity")

ggplot(adnc_app, aes(x = age, y = LVFFAEDVGSNK)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==1,], aes(x = age, y = LVFFAEDVGSNK)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_L", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==2,], aes(x = age, y = LVFFAEDVGSNK)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_M", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==3,], aes(x = age, y = LVFFAEDVGSNK)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_H", x = "age", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$age_group=='<60',], aes(x = adnc_fac, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_age_60", x = "adnc", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$age_group=='60-70',], aes(x = adnc_fac, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_age_60-70", x = "adnc", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$adnc==0,], aes(x = age_group, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_HC", x = "", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$adnc==1,], aes(x = age_group, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_L", x = "", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$adnc==2,], aes(x = age_group, y = LVFFAEDVGSNK)) +
  # geom_point(alpha = 0.5) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_M", x = "", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$adnc==3,], aes(x = age_group, y = LVFFAEDVGSNK)) +
  # geom_point(alpha = 0.5) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_H", x = "", y = "intensity")

# 控制时间
ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$age_group=='<60',], aes(x = adnc_fac, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_age_60", x = "adnc", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$age_group=='60-70',], aes(x = adnc_fac, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_age_60-70", x = "adnc", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$age_group=='70-80',], aes(x = adnc_fac, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_age_70-80", x = "adnc", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$age_group=='80-90',], aes(x = adnc_fac, y = LVFFAEDVGSNK)) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_age_80-90", x = "adnc", y = "intensity")

ggplot(peptide_APP_limma_t_im[peptide_APP_limma_t_im$age_group=='>90',], aes(x = adnc_fac, y = LVFFAEDVGSNK)) +
  # geom_point(alpha = 0.5) +
  geom_boxplot(alpha = 0.5) +  # 半透明 boxplot
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, color = "black") +  # 添加散点，并设置大小和颜色
  # ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) +
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_age_90", x = "adnc", y = "intensity")

# sex rin bank etal
# rin
ggplot(hc_app, aes(x = age, y = LVFFAEDVGSNK, color = rin)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  # ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_HC", x = "age", y = "intensity")

ggplot(adnc_app, aes(x = age, y = LVFFAEDVGSNK, color = rin)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==1,], aes(x = age, y = LVFFAEDVGSNK, color = rin)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_L", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==2,], aes(x = age, y = LVFFAEDVGSNK, color = rin)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_M", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==3,], aes(x = age, y = LVFFAEDVGSNK, color = rin)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_H", x = "age", y = "intensity")

# sex
ggplot(hc_app, aes(x = age, y = LVFFAEDVGSNK, color = as.factor(sex))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  # ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_HC", x = "age", y = "intensity")

ggplot(adnc_app, aes(x = age, y = LVFFAEDVGSNK, color = as.factor(sex))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==1,], aes(x = age, y = LVFFAEDVGSNK, color = as.factor(sex))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_L", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==2,], aes(x = age, y = LVFFAEDVGSNK, color = as.factor(sex))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_M", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==3,], aes(x = age, y = LVFFAEDVGSNK, color = as.factor(sex))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_H", x = "age", y = "intensity")

# bank
ggplot(hc_app, aes(x = age, y = LVFFAEDVGSNK, color = as.factor(bank))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) +
  # ylim(min(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE), max(peptide_APP_limma_t$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_HC", x = "age", y = "intensity")

ggplot(adnc_app, aes(x = age, y = LVFFAEDVGSNK, color = as.factor(bank))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==1,], aes(x = age, y = LVFFAEDVGSNK, color = as.factor(bank))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_L", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==2,], aes(x = age, y = LVFFAEDVGSNK, color = as.factor(bank))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_M", x = "age", y = "intensity")

ggplot(adnc_app[adnc_app$adnc==3,], aes(x = age, y = LVFFAEDVGSNK, color = as.factor(bank))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  ylim(min(adnc_app$LVFFAEDVGSNK, na.rm = TRUE), max(adnc_app$LVFFAEDVGSNK, na.rm = TRUE)) + 
  theme_minimal() +
  theme_bw() +
  labs(title = "LVFFAEDVGSNK_ADNC_H", x = "age", y = "intensity")


# 控制age
age_group <- c('<60','60-70','70-80','80-90','>90')

hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='<60'],
     main = 'APP_LVFFAEDVGSNK_age_60',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='60-70'],
     main = 'APP_LVFFAEDVGSNK_age_60-70',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='70-80'],
     main = 'APP_LVFFAEDVGSNK_age_70-80',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='70-80'&peptide_APP_limma_t_im$adnc==0],
     main = 'APP_LVFFAEDVGSNK_age_70-80_HC',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='70-80'&peptide_APP_limma_t_im$adnc==1],
     main = 'APP_LVFFAEDVGSNK_age_70-80_L',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='70-80'&peptide_APP_limma_t_im$adnc==2],
     main = 'APP_LVFFAEDVGSNK_age_70-80_M',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='70-80'&peptide_APP_limma_t_im$adnc==3],
     main = 'APP_LVFFAEDVGSNK_age_70-80_H',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='80-90'],
     main = 'APP_LVFFAEDVGSNK_age_80-90',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='80-90'&peptide_APP_limma_t_im$adnc==0],
     main = 'APP_LVFFAEDVGSNK_age_80-90_HC',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='80-90'&peptide_APP_limma_t_im$adnc==1],
     main = 'APP_LVFFAEDVGSNK_age_80-90_L',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='80-90'&peptide_APP_limma_t_im$adnc==2],
     main = 'APP_LVFFAEDVGSNK_age_80-90_M',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='80-90'&peptide_APP_limma_t_im$adnc==3],
     main = 'APP_LVFFAEDVGSNK_age_80-90_H',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='>90'],
     main = 'APP_LVFFAEDVGSNK_age_90',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='>90'&peptide_APP_limma_t_im$adnc==0],
     main = 'APP_LVFFAEDVGSNK_age_90_HC',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='>90'&peptide_APP_limma_t_im$adnc==1],
     main = 'APP_LVFFAEDVGSNK_age_90_L',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='>90'&peptide_APP_limma_t_im$adnc==2],
     main = 'APP_LVFFAEDVGSNK_age_90_M',
     xlab = '',
     ylab = '')
hist(peptide_APP_limma_t_im$LVFFAEDVGSNK[peptide_APP_limma_t_im$age_group=='>90'&peptide_APP_limma_t_im$adnc==3],
     main = 'APP_LVFFAEDVGSNK_age_90_H',
     xlab = '',
     ylab = '')



# # volcano plot
# # adnc
# peptide_all_adnc_limma_re$change = ifelse(peptide_all_adnc_limma_re$adj.P.Val < 0.05,  
#                                           ifelse(peptide_all_adnc_limma_re$logFC > 0 ,'Up','Down'),
#                                           'Stable')
# peptide_all_adnc_limma_re$label = ifelse(peptide_all_adnc_limma_re$adj.P.Val < 0.05,peptide_all_adnc_limma_re$gene_name,"")
# 
# ggplot(peptide_all_adnc_limma_re, aes(x = logFC, y = -log10(adj.P.Val), colour=change)) +
#   geom_point(alpha=0.4, size=2) +
#   scale_color_manual(values=c( "#546de5", "#d2dae2","#ff4757"))+
#   geom_vline(xintercept=0,lty=2,col="black",lwd=0.8) +
#   geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.8) +
#   labs(x="logFC",y="-log10(adj.P)")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid=element_line(color='white'),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12),
#         legend.position="right", 
#         legend.title = element_blank())+
#   geom_text_repel(data = peptide_all_adnc_limma_re,
#                   aes(x = logFC, y = -log10(adj.P.Val), label = label),
#                   size = 3, 
#                   segment.color = "black", 
#                   segment.size = 0.5,
#                   box.padding = 0.5,
#                   point.padding = 0.3,
#                   show.legend = FALSE)
# 
# 
# # braak
# peptide_all_braak_limma_re$change = ifelse(peptide_all_braak_limma_re$adj.P.Val < 0.05,  
#                                            ifelse(peptide_all_braak_limma_re$logFC > 0 ,'Up','Down'),
#                                            'Stable')
# peptide_all_braak_limma_re$label = ifelse(peptide_all_braak_limma_re$adj.P.Val < 0.05, peptide_all_braak_limma_re$gene_name,"")
# 
# ggplot(peptide_all_braak_limma_re, aes(x = logFC, y = -log10(adj.P.Val), colour=change)) +
#   geom_point(alpha=0.4, size=2) +
#   scale_color_manual(values=c( "#546de5", "#d2dae2","#ff4757"))+
#   geom_vline(xintercept=0,lty=2,col="black",lwd=0.8) +
#   geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.8) +
#   labs(x="logFC",y="-log10(adj.P)")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid=element_line(color='white'),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         legend.position="right", 
#         legend.title = element_blank())+
#   geom_text_repel(data = peptide_all_braak_limma_re, aes(x = logFC, 
#                                                          y = -log10(adj.P.Val), 
#                                                          label = label),
#                   size = 3,box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.8, "lines"), 
#                   segment.color = "black", 
#                   show.legend = FALSE)
# 
# 
# # abeta
# peptide_all_abeta_limma_re$change = ifelse(peptide_all_abeta_limma_re$adj.P.Val < 0.05,  
#                                            ifelse(peptide_all_abeta_limma_re$logFC > 0 ,'Up','Down'),
#                                            'Stable')
# peptide_all_abeta_limma_re$label = ifelse(peptide_all_abeta_limma_re$adj.P.Val < 0.05, peptide_all_abeta_limma_re$gene_name,"")
# 
# ggplot(peptide_all_abeta_limma_re, aes(x = logFC, y = -log10(adj.P.Val), colour=change)) +
#   geom_point(alpha=0.4, size=2) +
#   scale_color_manual(values=c( "#546de5", "#d2dae2","#ff4757"))+
#   geom_vline(xintercept=0,lty=2,col="black",lwd=0.8) +
#   geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.8) +
#   labs(x="logFC",y="-log10(adj.P)")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid=element_line(color='white'),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         legend.position="right", 
#         legend.title = element_blank())+
#   geom_text_repel(data = peptide_all_abeta_limma_re, aes(x = logFC, 
#                                                          y = -log10(adj.P.Val), 
#                                                          label = label),
#                   size = 3,box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.8, "lines"), 
#                   segment.color = "black", 
#                   show.legend = FALSE)
# 
# 
# # cscore
# peptide_all_cscore_limma_re$change = ifelse(peptide_all_cscore_limma_re$adj.P.Val < 0.05,  
#                                             ifelse(peptide_all_cscore_limma_re$logFC > 0 ,'Up','Down'),
#                                             'Stable')
# peptide_all_cscore_limma_re$label = ifelse(peptide_all_cscore_limma_re$adj.P.Val < 0.05, peptide_all_cscore_limma_re$gene_name,"")
# 
# ggplot(peptide_all_cscore_limma_re, aes(x = logFC, y = -log10(adj.P.Val), colour=change)) +
#   geom_point(alpha=0.4, size=2) +
#   scale_color_manual(values=c( "#546de5", "#d2dae2","#ff4757"))+
#   geom_vline(xintercept=0,lty=2,col="black",lwd=0.8) +
#   geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.8) +
#   labs(x="logFC",y="-log10(adj.P)")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust = 0.5), 
#         panel.grid=element_line(color='white'),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         legend.position="right", 
#         legend.title = element_blank())+
#   geom_text_repel(data = peptide_all_cscore_limma_re, aes(x = logFC, 
#                                                           y = -log10(adj.P.Val), 
#                                                           label = label),
#                   size = 3,box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.8, "lines"), 
#                   segment.color = "black", 
#                   show.legend = FALSE)
# 
