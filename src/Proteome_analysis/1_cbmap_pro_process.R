library(data.table)
library(ggplot2)
library(patchwork)
library(writexl)
library(dplyr)
library(tidyverse)
library(caret)
library(glmnet)


setwd('/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/pr_jingjie_966_20250214/XB07045B4DA_966samples_Preliminary_analysis_result/')


# data input
pro <- fread('MS_identified_information.txt')
sample_list <- as.data.frame(fread('/data/projects/sunly/sample.txt'))
pro_qu_raw <- as.data.frame(fread('Protein_Quant.tsv'))
sample_infor <- read.csv('/data/shared_data/China_Brain_MultiOmics/sample_information/final/CBMAP_sample_info_1187_final_20241101.csv')
sample_765 <- fread('processed_data/all_sample.txt')
batch_info <- data.frame(fread('processed_data/batch_info.txt'))
lbd <- as.data.frame(fread('/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/processed/CBMAP_sample_info_1187_final_20250528_LBD.csv'))


# sample information prepare
mix_sample <- paste0('MIX',1:21)
repeat_sample <- c('repeat1','repeat3','repeat4','repeat5','repeat7')
repeat_nsample <- c('PTB425','XYA20221024','XYA20191227','A2024CBB001','A2022CBB082')
repeat_table <- data.frame(repeat_list = repeat_sample,
                           repeat_sample = repeat_nsample,
                           pairs = paste0('pair',1:5))
repeat_table <- data.frame(
  pairs = rep(repeat_table$pairs, each = 2),
  value = c(rbind(repeat_table$repeat_list, repeat_table$repeat_sample))
)

sample_list[which(sample_list$ID == 'repeat2'),1:4] <- c('PTB087', 'PTB087','human','pumc')
sample_list <- sample_list[which(sample_list$ID %in% c(sample_765$jingjie_ID, mix_sample, repeat_sample)),]
sample_list$type_bank <- ifelse(sample_list$type_bank=='fetal', 'csu', sample_list$type_bank)
sample_infor$id <- ifelse(sample_infor$bank=='zju', paste0('A',sample_infor$id), sample_infor$id)
sample_list$sex <- sample_infor[match(sample_list$ID, sample_infor$id),'sex_male']
sample_list$ip_batch <- batch_info[match(sample_list$ID, batch_info$sample),'IP']
sample_list$enzy_batch <- batch_info[match(sample_list$ID, batch_info$sample),'enzymolysis']
sample_list$sample_type <- ifelse(sample_list$ID%in%repeat_nsample, 'original',
                                  ifelse(sample_list$ID%in%repeat_sample, 'repeat', 'other'))
sample_list$mix <- ifelse(sample_list$ID%in%mix_sample, 'MIX' ,'Non-MIX') 
sample_list$pair <- repeat_table$pairs[match(sample_list$ID, repeat_table$value)]
sample_list$pair[is.na(sample_list$pair)] <- 'other'

# protein data process
colnames(pro_qu_raw) <- sub("^[^_]+_([^\\.]+).*", "\\1", colnames(pro_qu_raw))
pro_qu_raw$PG.ProteinGroups <- sub(";.*", "", pro_qu_raw$PG.ProteinGroups)
colnames(pro_qu_raw)[colnames(pro_qu_raw) == 'repeat2'] <- 'PTB087'
colnames(pro_qu_raw)[1] <- 'protein'

pro_qu <- pro_qu_raw[,colnames(pro_qu_raw) %in% c('protein',sample_765$jingjie_ID,mix_sample,repeat_sample)]
pro_qu[,2:ncol(pro_qu)] <- lapply(pro_qu[,2:ncol(pro_qu)], function(x) {
  if (is.numeric(x)) {
    x[is.nan(x)] <- NA
  }
  return(x)
})
pro_qu$missing_rate <- apply(pro_qu[,2:ncol(pro_qu)], 1, function(x) mean(is.na(x)))
pro_qu <- pro_qu[,c(1,ncol(pro_qu),2:(ncol(pro_qu)-1))]


#################### PCA by using proteins without missing, 5258 proteins
pro_qu_select <- pro_qu
pro_qu_pca <- pro_qu_select[which(pro_qu_select$missing_rate == 0),] #5258

pca_mix <- function(data){
  rownames(data) <- data$protein
  data_t <- data.frame(t(data[,c(3:ncol(data))]))
  
  pca <- prcomp(data_t, scale. = TRUE)
  pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
  pca_df$type <- sample_list[match(rownames(pca_df),sample_list$ID), 'type_bank']
  pca_df$type <- as.factor(pca_df$type)
  pca_df$sex <- sample_list[match(rownames(pca_df),sample_list$ID), 'sex']
  pca_df$sex <- as.factor(pca_df$sex)
  pca_df$ip_batch <- sample_list[match(rownames(pca_df),sample_list$ID), 'ip_batch']
  pca_df$ip_batch <- factor(pca_df$ip_batch, levels = c(paste0('B',1:21)))
  pca_df$enzy_batch <- sample_list[match(rownames(pca_df),sample_list$ID), 'enzy_batch']
  pca_df$enzy_batch <- factor(pca_df$enzy_batch, levels = c(paste0('A',1:4)))
  pca_df$sample_type <- sample_list[match(rownames(pca_df),sample_list$ID), 'sample_type']
  pca_df$sample_type <- factor(pca_df$sample_type, levels = c('other','repeat','original'))
  pca_df$pair <- sample_list[match(rownames(pca_df),sample_list$ID), 'pair']
  pca_df$pair <- factor(pca_df$pair, levels = c(paste0('pair',1:5), 'other'))
  pca_df$Mix <- sample_list[match(rownames(pca_df),sample_list$ID), 'mix']
  pca_df_bank <- pca_df[pca_df$type%in%c('csu','pumc','zju'),]
  pca_df_bank$type <- factor(pca_df_bank$type, levels = c('csu','pumc','zju'))
  pca_df_sex <- pca_df[!is.na(pca_df$sex),]
  pca_df_sex$sex_male <- ifelse(pca_df_sex$sex==1,'male','female')

  
  pca_dist <- sqrt(pca$x[,1]^2 + pca$x[,2]^2)
  pca_outliers <- which(pca_dist > mean(pca_dist) + 3*sd(pca_dist))
  outlier_samples <- rownames(pca_df)[pca_outliers] 
  print(outlier_samples)
  
  common_theme <- theme_bw() +
    theme(
      panel.grid = element_line(color = "white"), 
      plot.title = element_text(size = 22, face = "bold"),
      legend.title = element_text(size = 14),  # 图例标题字体
      legend.text = element_text(size = 14),   # 图例内容字体
      axis.text.x = element_text(size = 14),   # x 轴刻度字体
      axis.text.y = element_text(size = 14),   # y 轴刻度字体
      axis.title.x = element_text(size = 14),  # x 轴标题字体
      axis.title.y = element_text(size = 14)   # y 轴标题字体
    )
  # repeat
  pca_repeat <- ggplot() +
    geom_point(data = subset(pca_df, sample_type == "other"),
               aes(x = PC1, y = PC2, shape = sample_type, color = pair),
               size = 2, alpha = 1) +
    geom_point(data = subset(pca_df, sample_type != "other"),
               aes(x = PC1, y = PC2, shape = sample_type, color = pair),
               size = 2, alpha = 1) +
    scale_color_manual(values = c(
      "pair1" = "#A6CEE3", "pair2" = "#B2DF8A", "pair3" = "#FDBF6F",
      "pair4" = "#CAB2D6", "pair5" = "#FF9999", "other" = "lightgray")) +
    stat_ellipse(data = pca_df, aes(x = PC1, y = PC2, color = pair), level = 0.95) +
    guides(
      shape = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    common_theme
  print(pca_repeat)
  
  # mix
  pca_mix <- ggplot() +
    geom_point(data = subset(pca_df, Mix == "Non-MIX"),
               aes(x = PC1, y = PC2, color = Mix),
               size = 2, alpha = 0.9) +
    geom_point(data = subset(pca_df, Mix == "MIX"),
               aes(x = PC1, y = PC2, color = Mix),
               size = 2, alpha = 0.9) +
    scale_color_manual(values = c("MIX" = "#FF9999", "Non-MIX" = "lightgray")) +
    stat_ellipse(data = pca_df, aes(x = PC1, y = PC2, color = Mix), level = 0.95) +
    labs(color = "Type") +
    common_theme
  print(pca_mix)
  
  # bank
  pca_type <- ggplot(pca_df_bank, aes(PC1, PC2, color = type)) +
    geom_point(size = 2) +
    stat_ellipse(level = 0.95) +
    labs(color = "Bank") +
    common_theme
  print(pca_type)
 
  # ip batch 
  pca_ip_batch <- ggplot(pca_df, aes(PC1, PC2, color = ip_batch)) +
    geom_point(size = 2) +
    stat_ellipse(level = 0.95) +
    common_theme
  print(pca_ip_batch)
  
  # enzyme batch
  pca_enzy_batch <- ggplot(pca_df, aes(PC1, PC2, color = enzy_batch)) +
    geom_point(size = 2) +
    stat_ellipse(level = 0.95) +
    common_theme
  
  # sex
  pca_sex <- ggplot(pca_df_sex, aes(PC1, PC2, color = sex_male)) +
    geom_point(size = 2) +
    stat_ellipse(level = 0.95) +
    common_theme
  print(pca_sex)
  
  pca <- pca_repeat + pca_mix + pca_sex + pca_type + pca_enzy_batch + pca_ip_batch + plot_layout(ncol = 2, nrow = 3) 
  print(pca)
  
  ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/sample_info/pca.pdf',
         plot = pca,
         width = 10,
         height = 10)
  
  return(list(pca_df = pca_df, 
              outlier_samples = outlier_samples,
              p_type = p_type, 
              p_ip_batch = p_ip_batch, 
              p_enzy_batch = p_enzy_batch,
              p_sex = p_sex,
              p_pca = pca))
}

result <- pca_mix(pro_qu_pca)
pca_df <- result$pca_df
out_sample <- result$outlier_samples
result$p_type
result$p_ip_batch
result$p_enzy_batch
result$p_sex

ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/sample_info/pca_sample.pdf', 
       plot = result$p_type, 
       width = 8, 
       height = 6)
ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/sample_info/pca_sex.pdf', 
       plot = result$p_sex, 
       width = 8, 
       height = 6)
ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/sample_info/ip_batch.pdf', 
       plot = result$p_ip_batch, 
       width = 8, 
       height = 6)
ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/sample_info/enzy_batch.pdf', 
       plot = result$p_enzy_batch, 
       width = 8, 
       height = 6)


# outlier sample remove
pro_qu_filter <- pro_qu_select[,-which(colnames(pro_qu_select)=='A2021CBB015')] #去除最偏的一个样本
keep_cols <- intersect(colnames(pro_qu_filter), sample_765$jingjie_ID[sample_765$CBMAP_profile_1187 != ''])
pro_qu_filter <- pro_qu_filter[,c('protein','missing_rate',keep_cols)]
pro_qu_filter$missing_rate <- apply(pro_qu_filter[,3:ncol(pro_qu_filter)], 1, function(x) mean(is.na(x)))
pro_qu_filter[,'gene_name'] <- pro[match(pro_qu_filter$protein, pro$`Protein accession`), 'Gene name']
pro_qu_filter <- pro_qu_filter[,c(ncol(pro_qu_filter),1,2,3:(ncol(pro_qu_filter)-1))]
write.table(pro_qu_filter,
            '/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/processed/processed_data_dedupe/human_brain_pro_729.txt', 
            sep = '\t', quote = F, row.names = F)


# protein missing rate 
missing_rate_df <- data.frame(missing_rate = pro_qu_filter$missing_rate)

p_missing <- ggplot(missing_rate_df, aes(x = missing_rate)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "#6BAED6", alpha = 0.7) +
  geom_text(stat = "bin", binwidth=0.05, aes(label = after_stat(count)), vjust = -0.5, size = 3) +
  labs(title = NULL,
       x = "Missing Rate",
       y = "Frequency") +
  theme_bw() + 
  theme(panel.grid = element_line(color = 'white'),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))
p_missing

ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/sample_info/protein_missing_rate.pdf', 
       plot = p_missing, 
       width = 8, 
       height = 8)


# log
pro_filter <- pro_qu_filter
pro_filter[,4:ncol(pro_filter)] <- log2(pro_filter[,4:ncol(pro_filter)]+1)


# remove proteins with missing rate > 50%
pro_filter <- pro_filter[which(pro_filter$missing_rate < 0.5),] 
rownames(pro_filter) <- paste0(pro_filter$gene_name,'_',pro_filter$protein)

write.table(pro_filter,
            '/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/processed/processed_data/human_brain_pro_50_log_729.txt',
            sep = '\t', row.names = F, quote = F)


# imputation and median normalization
library(impute)
set.seed(20250514)
pro_filter_im <- impute.knn(as.matrix(pro_filter[,4:ncol(pro_filter)]), k = 10, rng.seed = 2025)$data %>% as.data.frame()
pro_filter_im_norm <- sweep(pro_filter_im, 2, apply(pro_filter_im, 2, median, na.rm = TRUE), FUN = "-")


# rint
rint_transform <- function(x) {
  n <- length(x)
  ranks <- rank(x, ties.method = "average")
  qnorm((ranks - 0.5) / n)
}

pro_filter_im_norm_rint <- apply(pro_filter_im_norm, 1, rint_transform)
pro_filter_im_norm_rint <- t(pro_filter_im_norm_rint)

write.table(pro_filter_im_norm_rint, 
            '/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/processed/processed_data/to_szy/cbmap_pro_qc.txt', 
            sep = '\t', quote = F, row.names = T)

###################### sex match
pro_filter_im_norm_rint <- read.delim("/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/processed/processed_data/human_brain_pro_729_50_norm_rint.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
pro_50_raw <- data.frame(pro_filter_im_norm_rint)
pro_50_raw$gene_protein <- rownames(pro_50_raw)
pro_50_raw$gene_name <- sub("_.*", "", pro_50_raw$gene_protein)
pro_50_raw$protein <- sub(".*_", "", pro_50_raw$gene_protein)
pro_50_raw <- pro_50_raw[,c(730:732,1:729)]

sample_id <- intersect(sample_infor$id, colnames(pro_50_raw))
pro <- pro_50_raw[,c('protein',sample_id)]
rownames(pro) <- pro$protein
pro <- data.frame(t(pro[,-1]))

sample_sex <- sample_infor[sample_infor$id%in%sample_id,]

# training data and test data prepare
set.seed(20250825) 
trainIndex <- createDataPartition(sample_sex$sex_male, p = 0.7, list = FALSE)
sample_train <- sample_sex[trainIndex, ]
sample_test  <- sample_sex[-trainIndex, ]

pro_train <- pro[rownames(pro)%in%sample_train$id,]
pro_test <- pro[rownames(pro)%in%sample_test$id,]

x_train <- as.matrix(pro_train)
y_train <- as.factor(sample_train$sex_male) 
x_test  <- as.matrix(pro_test)
y_test  <- as.factor(sample_test$sex_male)

# LASSO logistic regression
cvfit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)
best_lambda <- cvfit$lambda.min
pred_prob <- predict(cvfit, s = best_lambda, newx = x_test, type = "response")
pred_class <- ifelse(pred_prob > 0.5, levels(y_train)[2], levels(y_train)[1])
accuracy <- mean(pred_class == y_test)
cat("LASSO Logistic 回归 准确率:", accuracy, "\n")

x_all <- as.matrix(pro)
y_all <- as.factor(sample_sex$sex_male)

pred_prob_all <- predict(cvfit, s = best_lambda, newx = x_all, type = "response")
pred_class_all <- ifelse(pred_prob_all > 0.5, levels(y_train)[2], levels(y_train)[1])

comparison <- data.frame(
  SampleID = rownames(pro), 
  TrueSex = as.character(y_all),
  PredSex = as.character(pred_class_all),
  Prob = as.numeric(pred_prob_all)
)

# sex mismatch sample
mismatch <- comparison[comparison$TrueSex != comparison$PredSex, ]

# hist plot
comparison <- comparison %>%
  mutate(mismatch = PredSex != TrueSex)

ggplot(comparison, aes(x = Prob)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black") +
  geom_text(data = mismatch, 
            aes(x = Prob, y = 5, label = SampleID), # y=5 可以调整到合适高度
            angle = 90, vjust = -0.5, hjust = 0, size = 3, color = "red") +
  geom_vline(xintercept = mismatch$Prob, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(x = "Predicted probability", y = "Count",
       title = "Histogram of predicted probabilities with mismatches")

#################### cor
calculate_correlation_matrix <- function(data) {
  
  pro_repeat <- data[,repeat_sample]
  pro_repeat_sample <- data[,repeat_nsample]
  pro_filter_all <- data[,-which(colnames(data) %in% c('protein','missing_rate','gene_name',repeat_sample,repeat_nsample, mix_sample))]
  
  # 计算repeat与对应样本相关性
  cor_repeat <- list()
  for (i in 1:ncol(pro_repeat)) { 
    x <- pro_repeat[[i]]  
    y <- pro_repeat_sample[[i]]
    
    lm_model <- lm(y ~ x)
    summary_lm <- summary(lm_model)
    
    if (summary_lm$coefficients[2, 4] < 0.05) {
      cor_method <- "pearson" 
    } else {
      cor_method <- "spearman"  
    }
    
    cor_repeat[[i]] <- data.frame(
      x = colnames(pro_repeat)[i],
      y = colnames(pro_repeat_sample)[i],
      cor = cor(x, y, use = "complete.obs", method = cor_method),
      method = cor_method
    )
  }
  cor_repeat <- do.call(rbind, cor_repeat)
  cor_repeat$group <- 'repeat-original'
  
  # 计算repeat与其他样本之间的相关性
  cor_repeat_all <- list()
  for (i in 1:ncol(pro_repeat)) {
    x <- pro_repeat[[i]]
    cor_repeat_other <- list()
    for (j in 1:ncol(pro_filter_all)) {
      y <- pro_filter_all[[j]]
      
      lm_model <- lm(y ~ x)
      summary_lm <- summary(lm_model)
      
      if (summary_lm$coefficients[2, 4] < 0.05) {
        cor_method <- "pearson" 
      } else {
        cor_method <- "spearman"  
      }
      
      cor_repeat_other[[j]] <- data.frame(
        x = colnames(pro_repeat)[i],
        y = colnames(pro_filter_all)[j],
        cor = cor(x, y, use = "complete.obs", method = cor_method),
        method = cor_method
      )
    }
    cor_repeat_all[[i]] <- do.call(rbind, cor_repeat_other)
  }
  cor_repeat_all <- do.call(rbind, cor_repeat_all)
  cor_repeat_all$group <- 'repeat-other'
  
  # 合并两个相关性矩阵
  cor_matrix <- rbind(cor_repeat, cor_repeat_all)
  # box plot
  cor_matrix_plot <- cor_matrix[which(cor_matrix$cor > 0.9),]
  p <- ggplot(cor_matrix_plot, aes(x = group, y = cor, fill = group)) +
    geom_boxplot() +
    theme_minimal() +
    theme_bw() +
    theme(panel.grid=element_line(color='white'),
          axis.text.x = element_text(size = 14),  # 设置x轴标签大小
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = 'none')
  print(p)
  return(list(cor_matrix = cor_matrix,
              boxplot = p)) 
}

pro_qu_log <- pro_qu
pro_qu_log[,c(3:ncol(pro_qu_log))] <- log(pro_qu_log[,c(3:ncol(pro_qu_log))] + 1)
cor_result_log <- calculate_correlation_matrix(pro_qu_log) #log(intensity+1)，缺失率小于50%蛋白
cor_matrix_log <- cor_result_log$cor_matrix

p <- p_missing / plot_spacer() / cor_result_log$boxplot +
  plot_layout(ncol = 1, heights = c(1, 0.15, 1))
p

ggsave('/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/plot/sample_info/missing_cor.pdf', 
       plot = p, 
       width = 10, 
       height = 10)

