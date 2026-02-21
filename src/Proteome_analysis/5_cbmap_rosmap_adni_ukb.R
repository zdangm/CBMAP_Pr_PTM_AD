library(impute)
library(writexl)
library(openxlsx)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpointdensity)
library(viridis)
library(MASS)
library(ggpubr)


source('/share/home/sunly/Rscript/cbmap_pro/cbmap_586/2_sample_process.R')
out_dir <- '/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/'

###################### cbmap limma
adnc_id <- intersect(colnames(pro_50_raw), sample_adnc$id)
pro_cbmap_limma <- pro_50_raw[,colnames(pro_50_raw)%in%adnc_id]
pro_cbmap_limma <- pro_cbmap_limma[, match(adnc_id, colnames(pro_cbmap_limma))]
sample_adnc <- sample_adnc[sample_adnc$id%in%adnc_id,]
sample_adnc <- sample_adnc[match(adnc_id, sample_adnc$id), ]

design <- model.matrix( ~ adnc_num + age + sex_male + PMD + RIN, data = sample_adnc)
fit <- lmFit(pro_cbmap_limma, design)
fit <- eBayes(fit)
adnc_num_lim_re <- topTable(fit, coef = "adnc_num", number = Inf, adjust.method = "BH")
adnc_num_lim_re$gene_protein <- rownames(adnc_num_lim_re)
adnc_num_lim_re$genename <- pro_50_raw$gene_name[match(adnc_num_lim_re$gene_protein, pro_50_raw$gene_protein)] 
adnc_num_lim_re$protein <- pro_50_raw$protein[match(adnc_num_lim_re$gene_protein, pro_50_raw$gene_protein)] 

##################### rosmap data process
# pro_rosmap<-read.csv("/data/shared_data/ROSMAP/Processed/proteomics/TMT_quantitation/Round1_400individuals/C2.median_polish_corrected_log2_abundanceRatioCenteredOnMedianOfBatchMediansPerProtein-8817x400.csv", header = TRUE, sep = ",",check.names = F)
# rownames(pro_rosmap) <- pro_rosmap[,1]
# pro_rosmap <- pro_rosmap[,-1]
# set.seed(20250514)
# pro_rosmap <- impute.knn(as.matrix(pro_rosmap), k=10, rng.seed= 2025)$data %>% as.data.frame()
# rint_transform <- function(x) {
#   n <- length(x)
#   ranks <- rank(x, ties.method = "average")
#   qnorm((ranks - 0.5) / n)
# }
# pro_rosmap <- apply(pro_rosmap, 1, rint_transform)
# pro_rosmap <- data.frame(t(pro_rosmap))
# 
# rosmap_adnc <- read.xlsx('/data/shared_data/ROSMAP/Processed/metadata/dataset_1495_cross-sectional_02-16-2025.xlsx')
# ROSMAP_assay_proteomics_info<-as.data.frame(fread("/data/shared_data/ROSMAP/Processed/metadata/ROSMAP_assay_proteomics_TMTquantitation_metadata.csv",check.names = F))
# ROSMAP_assay_proteomics_info <- ROSMAP_assay_proteomics_info[ROSMAP_assay_proteomics_info$controlType=='',]
# ROSMAP_assay_proteomics_info$individual_id <- sub(".*\\.", "", ROSMAP_assay_proteomics_info$specimenID)
# rin <- as.data.frame(fread('/data/shared_data/ROSMAP/Processed/metadata/ROSMAP_assay_rnaSeq_metadata.csv'))
# rna_id_match <- as.data.frame(fread('/data/shared_data/ROSMAP/Processed/metadata/ROSMAP_biospecimen_metadata.csv'))
# rin$individualID <- rna_id_match$individualID[match(rin$specimenID, rna_id_match$specimenID)]
# 
# sampleinfo_rosmap<-as.data.frame(fread("/data/shared_data/ROSMAP/Processed/metadata/ROSMAP_clinical.csv",check.names = F))
# sampleinfo_rosmap$sample_id <- ROSMAP_assay_proteomics_info[match(sampleinfo_rosmap$individualID,ROSMAP_assay_proteomics_info$individual_id),"batchChannel"]
# sampleinfo_rosmap$RIN <- rin$RIN[match(sampleinfo_rosmap$individualID,rin$individualID)]
# sampleinfo_rosmap$adnc <- rosmap_adnc[match(sampleinfo_rosmap$projid, rosmap_adnc$projid),'ADNC_4level']
# sampleinfo_rosmap$dementia <- ifelse(sampleinfo_rosmap$dcfdx_lv%in%c(4,5,6),1,0)
# # write.table(sampleinfo_rosmap,'processed_data/tmp/sample_info_rosmap.txt', sep = '\t', quote = F, row.names = F)
# 
# pro_rosmap_info <- data.frame(
#   sample = colnames(pro_rosmap)
# )
# pro_rosmap_info$age <- sampleinfo_rosmap[match(pro_rosmap_info$sample, sampleinfo_rosmap$sample_id),"age_death"]
# pro_rosmap_info$age<-ifelse(pro_rosmap_info$age=="90+","90",pro_rosmap_info$age)
# pro_rosmap_info$age<-as.numeric(pro_rosmap_info$age)
# pro_rosmap_info$sex <- sampleinfo_rosmap[match(pro_rosmap_info$sample, sampleinfo_rosmap$sample_id),"msex"]
# pro_rosmap_info$pmi <- sampleinfo_rosmap[match(pro_rosmap_info$sample, sampleinfo_rosmap$sample_id),"pmi"]
# pro_rosmap_info$rin <- sampleinfo_rosmap[match(pro_rosmap_info$sample, sampleinfo_rosmap$sample_id),"RIN"]
# pro_rosmap_info$adnc <- sampleinfo_rosmap[match(pro_rosmap_info$sample, sampleinfo_rosmap$sample_id),"adnc"]
# pro_rosmap_info$braak <- sampleinfo_rosmap[match(pro_rosmap_info$sample, sampleinfo_rosmap$sample_id),"braaksc"]
# pro_rosmap_info$dementia <- sampleinfo_rosmap[match(pro_rosmap_info$sample, sampleinfo_rosmap$sample_id),'dementia']
# 
# pro_rosmap_info <- pro_rosmap_info[!is.na(pro_rosmap_info$adnc),]
# set.seed(20250531)
# rownames(pro_rosmap_info) <- pro_rosmap_info$sample
# pro_rosmap_info <- impute.knn(as.matrix(pro_rosmap_info[,2:ncol(pro_rosmap_info)]), k=10, rng.seed = 20250531)$data %>% as.data.frame()
# pro_rosmap_info$sample <- rownames(pro_rosmap_info)
# 
# # rosmap limma
# rosmap_adnc_id <- intersect(pro_rosmap_info$sample, colnames(pro_rosmap))
# pro_rosmap <- pro_rosmap[,colnames(pro_rosmap)%in%pro_rosmap_info$sample]
# pro_rosmap <- pro_rosmap[,match(rosmap_adnc_id, colnames(pro_rosmap))]
# pro_rosmap_info <- pro_rosmap_info[match(rosmap_adnc_id, pro_rosmap_info$sample),]
# 
# # adnc
# design <- model.matrix( ~ adnc + age + sex + pmi + rin, data = pro_rosmap_info)
# fit <- lmFit(pro_rosmap, design)
# fit <- eBayes(fit)
# adnc_num_lim_rosmap <- topTable(fit, coef = "adnc", number = Inf, adjust.method = "BH")
# adnc_num_lim_rosmap$pro <- rownames(adnc_num_lim_rosmap)
# adnc_num_lim_rosmap$gene_name <- sub("\\|.*", "", adnc_num_lim_rosmap$pro)
# adnc_num_lim_rosmap$pro <- sub(".*\\|", "", adnc_num_lim_rosmap$pro)
# 
# write.table(adnc_num_lim_rosmap, paste0(out_dir, 'table/brain_blood_csf/adnc_cbmap_rosmap_limma.txt'), sep = '\t', quote = F, row.names = F)

ad_banner_rosmap <- read.xlsx('dea_rosmap_banner_advscontrol.xlsx')
ad_banner_rosmap <- read.xlsx('processed_data/tmp/ad_banner_rosmap_dea.xlsx')

# ukb limma
pro_ukb <- as.data.frame(fread('/data/projects/AD_PWAS_UKB/AD_epi/data/pro_epi.txt'))
pro_ukb$ID <- paste0('X', pro_ukb$ID)
rownames(pro_ukb) <- pro_ukb$ID
colnames(pro_ukb)[12:ncol(pro_ukb)] <- substring(colnames(pro_ukb)[12:ncol(pro_ukb)],5)
pro_ukb <- pro_ukb[!is.na(pro_ukb$education),]
pro_ukb_cov <- pro_ukb[,c(1,4,5,7,8,10,11)]
colnames(pro_ukb_cov)[7] <- 'AD'

pro_ukb_limma <- pro_ukb[,c(12:ncol(pro_ukb))]
set.seed(20250514)
pro_ukb_limma <- impute.knn(as.matrix(pro_ukb_limma), k=10, rng.seed = 2025)$data %>% as.data.frame()
pro_ukb_limma <- data.frame(t(pro_ukb_limma))

design <- model.matrix( ~ AD + age + sex + education, data = pro_ukb_cov)
fit <- lmFit(pro_ukb_limma, design)
fit <- eBayes(fit)
ad_lim_ukb <- topTable(fit, coef = "AD", number = Inf, adjust.method = "BH")
ad_lim_ukb$pro <- rownames(ad_lim_ukb)
write.table(ad_lim_ukb, paste0(out_dir, 'table/brain_blood_csf/adnc_cbmap_ukb_limma.txt'), sep = '\t', quote = F, row.names = F)
# ad_lim_ukb <- as.data.frame(fread('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/table/brain_blood_csf/adnc_cbmap_ukb_limma.txt'))

# csf limma
csf_pro_other <- read.xlsx('csf_clinicalad_vs_control.xlsx')
csf_pro_other$Entrez.gene.symbol <- sub("\\|.*", "", csf_pro_other$Entrez.gene.symbol)

# 散点图
# CBMAP-ROSMAP
logfc_cbmap_rosmap <- data.frame(gene = intersect(adnc_num_lim_re$genename, ad_banner_rosmap$genename))
logfc_cbmap_rosmap$logfc_cbmap <- adnc_num_lim_re$logFC[match(logfc_cbmap_rosmap$gene, adnc_num_lim_re$genename)]
logfc_cbmap_rosmap$logfc_rosmap <- ad_banner_rosmap$logfc[match(logfc_cbmap_rosmap$gene, ad_banner_rosmap$genename)]
logfc_cbmap_rosmap$dep <- ifelse(logfc_cbmap_rosmap$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$P.Value<0.05], 'yes', 'no')
logfc_cbmap_rosmap$dep_fdr <- ifelse(logfc_cbmap_rosmap$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], 'yes', 'no')
logfc_cbmap_rosmap <- logfc_cbmap_rosmap[order(logfc_cbmap_rosmap$logfc_cbmap),]

cor.test(logfc_cbmap_rosmap$logfc_cbmap, logfc_cbmap_rosmap$logfc_rosmap)
cor.test(logfc_cbmap_rosmap$logfc_cbmap[logfc_cbmap_rosmap$dep=='yes'], logfc_cbmap_rosmap$logfc_rosmap[logfc_cbmap_rosmap$dep=='yes'])
cor.test(logfc_cbmap_rosmap$logfc_cbmap[logfc_cbmap_rosmap$dep_fdr=='yes'], logfc_cbmap_rosmap$logfc_rosmap[logfc_cbmap_rosmap$dep_fdr=='yes'])
table(logfc_cbmap_rosmap$dep)
table(logfc_cbmap_rosmap$dep_fdr)


# # CBMAP-ROSMAP
# logfc_cbmap_rosmap <- data.frame(gene = intersect(adnc_num_lim_re$genename, adnc_num_lim_rosmap$gene_name))
# logfc_cbmap_rosmap$logfc_cbmap <- adnc_num_lim_re$logFC[match(logfc_cbmap_rosmap$gene, adnc_num_lim_re$genename)]
# logfc_cbmap_rosmap$logfc_rosmap <- adnc_num_lim_rosmap$logFC[match(logfc_cbmap_rosmap$gene, adnc_num_lim_rosmap$gene_name)]
# logfc_cbmap_rosmap$dep <- ifelse(logfc_cbmap_rosmap$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$P.Value<0.05], 'yes', 'no')
# logfc_cbmap_rosmap$dep_fdr <- ifelse(logfc_cbmap_rosmap$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], 'yes', 'no')
# logfc_cbmap_rosmap <- logfc_cbmap_rosmap[order(logfc_cbmap_rosmap$logfc_cbmap),]
# 
# cor.test(logfc_cbmap_rosmap$logfc_cbmap, logfc_cbmap_rosmap$logfc_rosmap)
# cor.test(logfc_cbmap_rosmap$logfc_cbmap[logfc_cbmap_rosmap$dep=='yes'], logfc_cbmap_rosmap$logfc_rosmap[logfc_cbmap_rosmap$dep=='yes'])
# cor.test(logfc_cbmap_rosmap$logfc_cbmap[logfc_cbmap_rosmap$dep_fdr=='yes'], logfc_cbmap_rosmap$logfc_rosmap[logfc_cbmap_rosmap$dep_fdr=='yes'])
# table(logfc_cbmap_rosmap$dep)
# table(logfc_cbmap_rosmap$dep_fdr)

# fdr
p_fdr <- ggplot(logfc_cbmap_rosmap[logfc_cbmap_rosmap$dep_fdr == 'yes', ], 
                aes(x = logfc_cbmap, y = logfc_rosmap)) +
  geom_point(color = "#6495ED", size = 1.5, alpha = 0.7) +
  # 添加回归线和95% CI
  geom_smooth(method = "lm", color = "#6495ED", se = TRUE, fill = "skyblue1", alpha = 0.2) +
  # 添加相关性文字
  stat_cor(method = "pearson", 
           label.x = min(logfc_cbmap_rosmap$logfc_cbmap),
           label.y = max(logfc_cbmap_rosmap$logfc_rosmap),
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           size = 5) +
  labs(title = "CBMAP_BANNER/ROSMAP_FDR", 
       x = "logFC_CBMAP", 
       y = "logFC_BANNER/ROSMAP") +
  theme_minimal() +
  theme_bw()
p_fdr
ggsave(paste0(out_dir,'plot/bbc/cbmap_banner&rosmap_fdr.pdf'),
       p_fdr,
       height = 6,
       width = 6)

#p
p_p <- ggplot(logfc_cbmap_rosmap[logfc_cbmap_rosmap$dep=='yes',], aes(x = logfc_cbmap, y = logfc_rosmap)) +
  geom_point(color = "#6495ED", size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#6495ED", se = TRUE, fill = "skyblue1", alpha = 0.2) +
  # 添加相关性文字
  stat_cor(method = "pearson", 
           label.x = min(logfc_cbmap_rosmap$logfc_cbmap),
           label.y = max(logfc_cbmap_rosmap$logfc_rosmap),
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           size = 5) +
  labs(title = "CBMAP_BANNER/ROSMAP_P", x = "logFC_CBMAP", y = "logFC_BANNER/ROSMAP") +
  theme_minimal() +
  theme_bw()
p_p
ggsave(paste0(out_dir,'plot/bbc/cbmap_banner&rosmap_p.pdf'),
       p_p,
       height = 6,
       width = 6)

#all
dense_region <- logfc_cbmap_rosmap %>%
  filter(abs(logfc_cbmap) < 0.17, abs(logfc_rosmap) < 0.2)

sparse_region <- logfc_cbmap_rosmap %>%
  filter(abs(logfc_cbmap) >= 0.17 | abs(logfc_rosmap) >= 0.2)

dens <- kde2d(dense_region$logfc_cbmap, dense_region$logfc_rosmap, n = 200)

get_density <- function(x, y, dens) {
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  # 边界修正
  ix[ix < 1] <- 1; ix[ix > length(dens$x)] <- length(dens$x)
  iy[iy < 1] <- 1; iy[iy > length(dens$y)] <- length(dens$y)
  dens$z[cbind(ix, iy)]
}

dense_region <- dense_region %>%
  mutate(density = get_density(logfc_cbmap, logfc_rosmap, dens))

min_log_density <- min(log1p(dense_region$density), na.rm = TRUE)

color_scale <- viridis(100, option = "inferno")

log_dens_vals <- log1p(dense_region$density)
min_idx <- which.min(log_dens_vals)
min_color_idx <- round((log_dens_vals[min_idx] - min(log_dens_vals)) /
                         (max(log_dens_vals) - min(log_dens_vals)) * 99) + 1
min_color <- color_scale[min_color_idx]

# 合并数据
sparse_region$density <- NA
all_data <- rbind(sparse_region, dense_region)

# 计算相关系数
cor_test <- cor.test(all_data$logfc_cbmap, all_data$logfc_rosmap)
r_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 3)
cor_label <- paste0("R = ", r_val, ", p = ", p_val)

# 绘图
p_all <- ggplot() +
  geom_point(data = sparse_region,
             aes(x = logfc_cbmap, y = logfc_rosmap),
             color = min_color, alpha = 0.7, size = 1.5) +
  geom_point(data = dense_region,
             aes(x = logfc_cbmap, y = logfc_rosmap, color = log1p(density)),
             size = 1.5, alpha = 0.7) +
  scale_color_viridis(name = "log-density", option = "inferno") +
  geom_smooth(data = all_data,
              aes(x = logfc_cbmap, y = logfc_rosmap),
              method = "lm", se = TRUE,
              color = "#6495ED", fill = "skyblue1") +
  annotate("text",
           x = min(all_data$logfc_cbmap),
           y = max(all_data$logfc_rosmap),
           label = cor_label,
           hjust = 0, vjust = 1,
           size = 5, color = "black") +
  labs(x = "logFC_CBMAP", y = "logFC_BANNER/ROSMAP", title = "CBMAP_BANNER/ROSMAP_ALL") +
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

p_all
ggsave(paste0(out_dir,'plot/bbc/cbmap_banner&rosmap_all.pdf'),
       p_all,
       height = 6,
       width = 7)


# CBMAP-UKB
logfc_cbmap_ukb <- data.frame(gene = intersect(adnc_num_lim_re$genename, ad_lim_ukb$pro))
logfc_cbmap_ukb$logfc_cbmap <- adnc_num_lim_re$logFC[match(logfc_cbmap_ukb$gene, adnc_num_lim_re$genename)]
logfc_cbmap_ukb$logfc_ukb <- ad_lim_ukb$logFC[match(logfc_cbmap_ukb$gene, ad_lim_ukb$pro)]
logfc_cbmap_ukb$dep <- ifelse(logfc_cbmap_ukb$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$P.Value<0.05], 'yes', 'no')
logfc_cbmap_ukb$dep_fdr <- ifelse(logfc_cbmap_ukb$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], 'yes', 'no')
logfc_cbmap_ukb <- logfc_cbmap_ukb[order(logfc_cbmap_ukb$logfc_cbmap),]

cor.test(logfc_cbmap_ukb$logfc_cbmap, logfc_cbmap_ukb$logfc_ukb)
cor.test(logfc_cbmap_ukb$logfc_cbmap[logfc_cbmap_ukb$dep=='yes'], logfc_cbmap_ukb$logfc_ukb[logfc_cbmap_ukb$dep=='yes'])
cor.test(logfc_cbmap_ukb$logfc_cbmap[logfc_cbmap_ukb$dep_fdr=='yes'], logfc_cbmap_ukb$logfc_ukb[logfc_cbmap_ukb$dep_fdr=='yes'])
table(logfc_cbmap_ukb$dep_fdr)
table(logfc_cbmap_ukb$dep)
table(pro_ukb_cov$AD)

# fdr
cu_fdr <- ggplot(logfc_cbmap_ukb[logfc_cbmap_ukb$dep_fdr=='yes',], aes(x = logfc_cbmap, y = logfc_ukb)) +
  geom_point(color = "#6495ED", size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#6495ED", se = TRUE, fill = "skyblue1", alpha = 0.2) +
  # 添加相关性文字
  stat_cor(method = "pearson",
           label.x = min(logfc_cbmap_ukb$logfc_cbmap),
           label.y = max(logfc_cbmap_ukb$logfc_ukb),
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           size = 5) +
  labs(title = "CBMAP_UKB_FDR", x = "logFC_CBMAP", y = "logFC_UKB") +
  theme_minimal() +
  theme_bw()
cu_fdr

ggsave(paste0(out_dir,'plot/bbc/cbmap_ukb_fdr.pdf'),
       cu_fdr,
       height = 6,
       width = 6)


#p
cu_p <- ggplot(logfc_cbmap_ukb[logfc_cbmap_ukb$dep=='yes',], aes(x = logfc_cbmap, y = logfc_ukb)) +
  geom_point(color = "#6495ED", size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#6495ED", se = TRUE, fill = "skyblue1", alpha = 0.2) +
  # 添加相关性文字
  stat_cor(method = "pearson",
           label.x = min(logfc_cbmap_ukb$logfc_cbmap),
           label.y = max(logfc_cbmap_ukb$logfc_ukb),
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           size = 5) +
  labs(title = "CBMAP_UKB_P", x = "logFC_CBMAP", y = "logFC_UKB") +
  theme_minimal() +
  theme_bw()
cu_p

ggsave(paste0(out_dir,'plot/bbc/cbmap_ukb_p.pdf'),
       cu_p,
       height = 6,
       width = 6)

# all
dense_region <- logfc_cbmap_ukb %>%
  filter(abs(logfc_cbmap) < 0.15, abs(logfc_ukb) < 0.12)

sparse_region <- logfc_cbmap_ukb %>%
  filter(abs(logfc_cbmap) >= 0.15 | abs(logfc_ukb) >= 0.12)

dens <- kde2d(dense_region$logfc_cbmap, dense_region$logfc_ukb, n = 200)

get_density <- function(x, y, dens) {
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  # 边界修正
  ix[ix < 1] <- 1; ix[ix > length(dens$x)] <- length(dens$x)
  iy[iy < 1] <- 1; iy[iy > length(dens$y)] <- length(dens$y)
  dens$z[cbind(ix, iy)]
}

dense_region <- dense_region %>%
  mutate(density = get_density(logfc_cbmap, logfc_ukb, dens))

min_log_density <- min(log1p(dense_region$density), na.rm = TRUE)

color_scale <- viridis(100, option = "inferno")

log_dens_vals <- log1p(dense_region$density)
min_idx <- which.min(log_dens_vals)
min_color_idx <- round((log_dens_vals[min_idx] - min(log_dens_vals)) /
                         (max(log_dens_vals) - min(log_dens_vals)) * 99) + 1
min_color <- color_scale[min_color_idx]

# 合并dense spare
sparse_region$density <- NA
all_data <- rbind(sparse_region, dense_region)

# 相关性
cor_test <- cor.test(all_data$logfc_cbmap, all_data$logfc_ukb)
r_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 3)
cor_label <- paste0("R = ", r_val, ", p = ", p_val)

# plot
cu_all <- ggplot() +
  geom_point(data = sparse_region,
             aes(x = logfc_cbmap, y = logfc_ukb),
             color = min_color, alpha = 0.7, size = 1.5) +
  geom_point(data = dense_region,
             aes(x = logfc_cbmap, y = logfc_ukb, color = log1p(density)),
             size = 1.5, alpha = 0.7) +
  geom_smooth(data = all_data,
              aes(x = logfc_cbmap, y = logfc_ukb),
              method = "lm", se = TRUE,
              color = "#6495ED", fill = "skyblue1") +
  annotate("text",
           x = min(all_data$logfc_cbmap),
           y = max(all_data$logfc_ukb),
           label = cor_label,
           hjust = 0, vjust = 1,
           size = 5, color = "black") +
  scale_color_viridis(name = "log-density", option = "inferno") +
  labs(x = "logFC_CBMAP", y = "logFC_UKB", title = 'CBMAP_UKB_ALL') +
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )
cu_all

ggsave(paste0(out_dir,'plot/bbc/cbmap_ukb_all.pdf'),
       cu_all,
       height = 6,
       width = 7)


# CBMAP-CSF
logfc_cbmap_csf_other <- data.frame(gene = intersect(csf_pro_other$Entrez.gene.symbol, adnc_num_lim_re$genename))
logfc_cbmap_csf_other$logfc_cbmap <- adnc_num_lim_re$logFC[match(logfc_cbmap_csf_other$gene, adnc_num_lim_re$genename)]
logfc_cbmap_csf_other$logfc_csf <- csf_pro_other$Beta[match(logfc_cbmap_csf_other$gene, csf_pro_other$Entrez.gene.symbol)]
logfc_cbmap_csf_other$dep <- ifelse(logfc_cbmap_csf_other$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$P.Value<0.05], 'yes', 'no')
logfc_cbmap_csf_other$dep_fdr <- ifelse(logfc_cbmap_csf_other$gene%in%adnc_num_lim_re$genename[adnc_num_lim_re$adj.P.Val<0.05], 'yes', 'no')
logfc_cbmap_csf_other <- logfc_cbmap_csf_other[order(logfc_cbmap_csf_other$logfc_cbmap),]

cor.test(logfc_cbmap_csf_other$logfc_cbmap, logfc_cbmap_csf_other$logfc_csf)
cor.test(logfc_cbmap_csf_other$logfc_cbmap[logfc_cbmap_csf_other$dep=='yes'], logfc_cbmap_csf_other$logfc_csf[logfc_cbmap_csf_other$dep=='yes'])
cor.test(logfc_cbmap_csf_other$logfc_cbmap[logfc_cbmap_csf_other$dep_fdr=='yes'], logfc_cbmap_csf_other$logfc_csf[logfc_cbmap_csf_other$dep_fdr=='yes'])
table(logfc_cbmap_csf_other$dep)
table(logfc_cbmap_csf_other$dep_fdr)

# fdr
cc_fdr <- ggplot(logfc_cbmap_csf_other[logfc_cbmap_csf_other$dep_fdr=='yes',], aes(x = logfc_cbmap, y = logfc_csf)) +
  geom_point(color = "#6495ED", size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#6495ED", se = TRUE, fill = "skyblue1", alpha = 0.2) +
  # 添加相关性文字
  stat_cor(method = "pearson",
           label.x = min(logfc_cbmap_csf_other$logfc_cbmap),
           label.y = max(logfc_cbmap_csf_other$logfc_csf),
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           size = 5) +
  labs(title = "CBMAP_ADNI_FDR", x = "logFC_CBMAP", y = "logFC_CSF") +
  theme_minimal() +
  theme_bw()
cc_fdr

ggsave(paste0(out_dir,'plot/bbc/cbmap_csf_fdr.pdf'),
       cc_fdr,
       height = 6,
       width = 6,
       dpi = 300)

#p
cc_p <- ggplot(logfc_cbmap_csf_other[logfc_cbmap_csf_other$dep=='yes',], aes(x = logfc_cbmap, y = logfc_csf)) +
  geom_point(color = "#6495ED", size = 1.5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "#6495ED", se = TRUE, fill = "skyblue1", alpha = 0.2) +
  # 添加相关性文字
  stat_cor(method = "pearson",
           label.x = min(logfc_cbmap_csf_other$logfc_cbmap),
           label.y = max(logfc_cbmap_csf_other$logfc_csf),
           aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
           size = 5) +
  labs(title = "CBMAP_ADNI_P", x = "logFC_CBMAP", y = "logFC_CSF") +
  theme_minimal() +
  theme_bw()
cc_p

ggsave(paste0(out_dir,'plot/bbc/cbmap_csf_p.pdf'),
       cc_p,
       height = 6,
       width = 6)

# all
dense_region <- logfc_cbmap_csf_other %>%
  filter(abs(logfc_cbmap) < 0.17, abs(logfc_csf) < 0.5)

sparse_region <- logfc_cbmap_csf_other %>%
  filter(abs(logfc_cbmap) >= 0.17 | abs(logfc_csf) >= 0.5)

dens <- kde2d(dense_region$logfc_cbmap, dense_region$logfc_csf, n = 200)

get_density <- function(x, y, dens) {
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  # 边界修正
  ix[ix < 1] <- 1; ix[ix > length(dens$x)] <- length(dens$x)
  iy[iy < 1] <- 1; iy[iy > length(dens$y)] <- length(dens$y)
  dens$z[cbind(ix, iy)]
}

dense_region <- dense_region %>%
  mutate(density = get_density(logfc_cbmap, logfc_csf, dens))

min_log_density <- min(log1p(dense_region$density), na.rm = TRUE)

color_scale <- viridis(100, option = "inferno")

log_dens_vals <- log1p(dense_region$density)
min_idx <- which.min(log_dens_vals)
min_color_idx <- round((log_dens_vals[min_idx] - min(log_dens_vals)) /
                         (max(log_dens_vals) - min(log_dens_vals)) * 99) + 1
min_color <- color_scale[min_color_idx]

# 合并dense spare
sparse_region$density <- NA
all_data <- rbind(sparse_region, dense_region)

# 相关性
cor_test <- cor.test(all_data$logfc_cbmap, all_data$logfc_csf)
r_val <- round(cor_test$estimate, 2)
p_val <- signif(cor_test$p.value, 3)
cor_label <- paste0("R = ", r_val, ", p = ", p_val)

cc_all <- ggplot() +
  geom_point(data = sparse_region,
             aes(x = logfc_cbmap, y = logfc_csf),
             color = min_color, alpha = 0.7, size = 1.5) +
  geom_point(data = dense_region,
             aes(x = logfc_cbmap, y = logfc_csf, color = log1p(density)),
             size = 1.5, alpha = 0.7) +
  geom_smooth(data = all_data,
              aes(x = logfc_cbmap, y = logfc_csf),
              method = "lm", se = TRUE,
              color = "#6495ED", fill = "skyblue1") +
  annotate("text",
           x = min(all_data$logfc_cbmap),
           y = max(all_data$logfc_csf),
           label = cor_label,
           hjust = 0, vjust = 1,
           size = 5, color = "black") +
  scale_color_viridis(name = "log-density", option = "inferno") +
  labs(x = "logFC_CBMAP", y = "logFC_CSF", title = 'CBMAP_ADNI_ALL') +
  theme_minimal() +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )
cc_all

ggsave(paste0(out_dir,'plot/bbc/cbmap_csf_all.pdf'),
       cc_all,
       height = 6,
       width = 7)

library(patchwork)

p <- p_fdr + p_p + p_all +
  plot_spacer() + plot_spacer() + plot_spacer() +  # 空白行
  cc_fdr + cc_p + cc_all +
  plot_spacer() + plot_spacer() + plot_spacer() +  # 空白行
  cu_fdr + cu_p + cu_all +
  plot_layout(ncol = 3, nrow = 5, heights = c(1, 0.1, 1, 0.1, 1))  # 注意行数变为5
p

ggsave(paste0(out_dir,'plot/bbc/all.pdf'),
       p,
       height = 18,
       width = 15)
ggsave(paste0(out_dir,'plot/bbc/all.png'),
       p,
       height = 18,
       width = 15,
       dpi = 300)

colnames(adnc_num_lim_re) <- paste0(colnames(adnc_num_lim_re),'.CBMAP')
colnames(ad_banner_rosmap) <- paste0(colnames(ad_banner_rosmap),'.ROSMAP')
colnames(ad_lim_ukb) <- paste0(colnames(ad_lim_ukb),'.UKB')
colnames(csf_pro_other)[8:9] <- c('P.Value','adj.P.Val')
colnames(csf_pro_other) <- paste0(colnames(csf_pro_other),'.ADNI')

adnc_cbmap_rosmap <- data.frame(genename = intersect(adnc_num_lim_re$genename.CBMAP, ad_banner_rosmap$genename.ROSMAP))
adnc_cbmap_rosmap$logfc.cbmap <- adnc_num_lim_re$logFC.CBMAP[match(adnc_cbmap_rosmap$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_rosmap$P.Value.cbmap <- adnc_num_lim_re$P.Value.CBMAP[match(adnc_cbmap_rosmap$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_rosmap$adj.P.Val.cbmap <- adnc_num_lim_re$adj.P.Val.CBMAP[match(adnc_cbmap_rosmap$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_rosmap$logfc.rosmap <- ad_banner_rosmap$logfc.ROSMAP[match(adnc_cbmap_rosmap$genename, ad_banner_rosmap$genename.ROSMAP)]
# adnc_cbmap_rosmap$P.Value.rosmap <- ad_banner_rosmap$P.Value.CBMAP[match(adnc_cbmap_rosmap$genename, ad_banner_rosmap$genename.ROSMAP)]
adnc_cbmap_rosmap$adj.P.Val.rosmap <- ad_banner_rosmap$adj.Pvalue.ROSMAP[match(adnc_cbmap_rosmap$genename, ad_banner_rosmap$genename.ROSMAP)]
adnc_cbmap_rosmap$dep_fdr <- ifelse(adnc_cbmap_rosmap$genename%in%adnc_num_lim_re$genename.CBMAP[adnc_num_lim_re$adj.P.Val.CBMAP<0.05], 'yes', 'no')
adnc_cbmap_rosmap$dep <- ifelse(adnc_cbmap_rosmap$genename%in%adnc_num_lim_re$genename.CBMAP[adnc_num_lim_re$P.Value.CBMAP<0.05], 'yes', 'no')


adnc_cbmap_ukb <- data.frame(genename = intersect(adnc_num_lim_re$genename.CBMAP, ad_lim_ukb$pro.UKB))
adnc_cbmap_ukb$logfc.cbmap <- adnc_num_lim_re$logFC.CBMAP[match(adnc_cbmap_ukb$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_ukb$P.Value.cbmap <- adnc_num_lim_re$P.Value.CBMAP[match(adnc_cbmap_ukb$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_ukb$adj.P.Val.cbmap <- adnc_num_lim_re$adj.P.Val.CBMAP[match(adnc_cbmap_ukb$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_ukb$logfc.ukb <- ad_lim_ukb$logFC.UKB[match(adnc_cbmap_ukb$genename, ad_lim_ukb$pro.UKB)]
adnc_cbmap_ukb$P.Value.ukb <- ad_lim_ukb$P.Value.UKB[match(adnc_cbmap_ukb$genename, ad_lim_ukb$pro.UKB)]
adnc_cbmap_ukb$adj.P.Val.ukb <- ad_lim_ukb$adj.P.Val.UKB[match(adnc_cbmap_ukb$genename, ad_lim_ukb$pro.UKB)]
adnc_cbmap_ukb$dep_fdr <- ifelse(adnc_cbmap_ukb$genename%in%adnc_num_lim_re$genename.CBMAP[adnc_num_lim_re$adj.P.Val.CBMAP<0.05], 'yes', 'no')
adnc_cbmap_ukb$dep <- ifelse(adnc_cbmap_ukb$genename%in%adnc_num_lim_re$genename.CBMAP[adnc_num_lim_re$P.Value.CBMAP<0.05], 'yes', 'no')


adnc_cbmap_adni <- data.frame(genename = intersect(adnc_num_lim_re$genename.CBMAP, csf_pro_other$Entrez.gene.symbol.ADNI))
adnc_cbmap_adni$logfc.cbmap <- adnc_num_lim_re$logFC.CBMAP[match(adnc_cbmap_adni$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_adni$P.Value.cbmap <- adnc_num_lim_re$P.Value.CBMAP[match(adnc_cbmap_adni$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_adni$adj.P.Val.cbmap <- adnc_num_lim_re$adj.P.Val.CBMAP[match(adnc_cbmap_adni$genename, adnc_num_lim_re$genename.CBMAP)]
adnc_cbmap_adni$logfc.ukb <- csf_pro_other$Beta.ADNI[match(adnc_cbmap_adni$genename, csf_pro_other$Entrez.gene.symbol.ADNI)]
adnc_cbmap_adni$P.Value.ukb <- csf_pro_other$P.Value.ADNI[match(adnc_cbmap_adni$genename, csf_pro_other$Entrez.gene.symbol.ADNI)]
adnc_cbmap_adni$adj.P.Val.ukb <- csf_pro_other$adj.P.Val.ADNI[match(adnc_cbmap_adni$genename, csf_pro_other$Entrez.gene.symbol.ADNI)]
adnc_cbmap_adni$dep_fdr <- ifelse(adnc_cbmap_adni$genename%in%adnc_num_lim_re$genename.CBMAP[adnc_num_lim_re$adj.P.Val.CBMAP<0.05], 'yes', 'no')
adnc_cbmap_adni$dep <- ifelse(adnc_cbmap_adni$genename%in%adnc_num_lim_re$genename.CBMAP[adnc_num_lim_re$P.Value.CBMAP<0.05], 'yes', 'no')


CBMAP_OTHER <- list(
  adnc_cbmap_rosmap = adnc_cbmap_rosmap,
  adnc_cbmap_ukb = adnc_cbmap_ukb,
  adnc_cbmap_ADNI = adnc_cbmap_adni
)

write_xlsx(CBMAP_OTHER, paste0(out_dir, 'table/brain_blood_csf/cbmap_bbc.xlsx'))


