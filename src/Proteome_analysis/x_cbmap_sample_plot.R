library(ggrepel)
library(ggplot2)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(impute)
library(UpSetR)
library(sf)


source('/share/home/sunly/Rscript/cbmap_pro/cbmap_586/2_sample_process.R')
out_dir <- '/data/projects/China_Brain_MultiOmics/humanBrain_protein/cbmap_pro_trait_asso/result/cbmap_586/'

# adnc braak age bank sex条形图
# adnc
sample_adnc <- sample_adnc[sample_adnc$id%in%colnames(pro_50_raw)[4:ncol(pro_50_raw)],]
adnc_counts <- as.data.frame(table(sample_adnc$adnc_num))
adnc_counts$adnc <- ifelse(adnc_counts$Var1==0,'HC',
                           ifelse(adnc_counts$Var1==1,'ADNC-L',
                                  ifelse(adnc_counts$Var1==2,'ADNC-M','ADNC-H')))
adnc_counts$adnc <- factor(adnc_counts$adnc, levels = c('HC','ADNC-L','ADNC-M','ADNC-H'))
my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(nrow(adnc_counts)) 
p1 <- ggplot(adnc_counts, aes(x = adnc, y = Freq, fill = adnc)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +  # 关闭图例
  geom_text(aes(label = Freq), vjust = -0.5, size = 5) +  # 标注数量
  scale_fill_manual(values = my_colors) +  # 使用自定义颜色
  theme_minimal(base_size = 14) +
  labs(title = "Sample information of ADNC", x = NULL, y = "Count") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_line(color = 'white'),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20)
  )
p1
ggsave(paste0(out_dir, 'plot/sample_info/adnc_sample_information.pdf'), 
       plot = p1, 
       width = 8, 
       height = 6)

# braak
sample_braak <- sample_braak[sample_braak$id%in%colnames(pro_50_raw)[4:ncol(pro_50_raw)],]
braak_counts <- table(sample_braak$braak_fac)
braak_counts <- as.data.frame(braak_counts)
colnames(braak_counts)[1] <- 'braak' 
my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(nrow(braak_counts)) 
p2 <- ggplot(braak_counts, aes(x = braak, y = Freq, fill = braak)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +  # 关闭图例
  geom_text(aes(label = Freq), vjust = -0.5, size = 5) +  # 标注数量
  scale_fill_manual(values = my_colors) +  # 使用自定义颜色
  theme_minimal(base_size = 14) +
  labs(title = "Sample information of Braak", x = NULL, y = "Count") +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid = element_line(color = 'white'),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20)
  )
p2
ggsave(paste0(out_dir, 'plot/sample_info/braak_sample_information.pdf'), 
       plot = p2, 
       width = 8, 
       height = 6)

# abeta
sample_abeta <- sample_abeta[sample_abeta$id%in%colnames(pro_50_raw)[4:ncol(pro_50_raw)],]
abeta_counts <- table(sample_abeta$A_beta_0_3)
abeta_counts <- as.data.frame(abeta_counts)
colnames(abeta_counts)[1] <- 'abeta' 
my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(nrow(abeta_counts)) 
p3 <- ggplot(abeta_counts, aes(x = abeta, y = Freq, fill = abeta)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +  # 关闭图例
  geom_text(aes(label = Freq), vjust = -0.5, size = 5) +  # 标注数量
  scale_fill_manual(values = my_colors) +  # 使用自定义颜色
  theme_minimal(base_size = 14) +
  labs(title = "Sample information of Abeta", x = NULL, y = "Count") +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid = element_line(color = 'white'),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20)
  )
p3
ggsave(paste0(out_dir,'plot/sample_info/abeta_sample_information.pdf'), 
       plot = p3, 
       width = 8, 
       height = 6)

# cscore
sample_cscore <- sample_cscore[sample_cscore$id%in%colnames(pro_50_raw)[4:ncol(pro_50_raw)],]
cscore_counts <- table(sample_cscore$C_cerad_0_3)
cscore_counts <- as.data.frame(cscore_counts)
colnames(cscore_counts)[1] <- 'cscore' 
my_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(nrow(cscore_counts)) 
p4 <- ggplot(cscore_counts, aes(x = cscore, y = Freq, fill = cscore)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +  # 关闭图例
  geom_text(aes(label = Freq), vjust = -0.5, size = 5) +  # 标注数量
  scale_fill_manual(values = my_colors) +  # 使用自定义颜色
  theme_minimal(base_size = 14) +
  labs(title = "Sample information of cerad", x = NULL, y = "Count") +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.grid = element_line(color = 'white'),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    plot.title = element_text(size = 20)
  )
p4
ggsave(paste0(out_dir,'plot/sample_info/cerad_sample_information.pdf'), 
       plot = p4, 
       width = 8, 
       height = 6)

# age
p5 <- ggplot(sample_adnc, aes(x = age, fill = after_stat(x))) +
  geom_histogram(binwidth = 10, color = "black") +
  scale_fill_gradient(low = "#C6DBEF", high = "#08306B") +
  geom_text(stat='bin', binwidth=10, aes(label=after_stat(count)), vjust=-0.5, size=5) +
  labs(x = "Age", y = "Count", fill = "Age", title = 'Age distribution') +
  theme_bw() + 
  theme(
    panel.grid = element_line(color = 'white'),
    legend.position = 'none',
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    plot.title = element_text(size = 20)
  )
p5
ggsave(paste0(out_dir,'plot/sample_info/age_sample_information.pdf'), 
       plot = p5, 
       width = 8, 
       height = 6)

# bank
sample_adnc$BANK <- ifelse(sample_adnc$bank=='zju','ZJU',
                           ifelse(sample_adnc$bank=='csu','CSU','PUMC'))
bank_counts <- sample_adnc %>%
  count(BANK) %>%
  mutate(prop = n / sum(n),
         label = paste0(BANK, "\n", n, " (", scales::percent(prop), ")")
         # label = paste0(BANK, ": (", n, ")")
  )
p6 <- ggplot(bank_counts, aes(x = "", y = prop, fill = BANK)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  labs(title = "") +
  theme_void() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_brewer(palette = "Set3") +
  theme(legend.position = "none",
        text = element_text(size = 8, family = "Arial") )
p6
ggsave(paste0(out_dir,'plot/sample_info/bank_sample_information.pdf'), 
       plot = p6, 
       width = 4, 
       height = 4)

# bank sex
sample_adnc$sex <- ifelse(sample_adnc$sex_male==1, 'Male', 'Female')
p7 <- ggplot(sample_adnc, aes(x = BANK, fill = sex)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", 
            aes(label = after_stat(count)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 4)  +  # 标注数量
  labs(title = NULL,
       x = NULL, y = "Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_bw() +
  theme(
    panel.grid = element_line(color = 'white'),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "none"
  )
p7
ggsave(paste0(out_dir,'plot/sample_info/bank_sex_sample_information.pdf'), 
       plot = p7, 
       width = 7, 
       height = 8)

# bank sex adnc
sample_adnc$disease <- ifelse(sample_adnc$adnc_num==0,'HC',
                              ifelse(sample_adnc$adnc_num==1,'ADNC-L',
                                     ifelse(sample_adnc$adnc_num==2,'ADNC-M','ADNC-H')))
sample_adnc$disease <- factor(sample_adnc$disease, levels = c('HC','ADNC-L','ADNC-M','ADNC-H'))
p8 <- ggplot(sample_adnc, aes(x = BANK, fill = sex)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", 
            aes(label = after_stat(count)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 4) +
  facet_wrap(~ disease) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  scale_fill_brewer(palette = "Pastel1") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14)
  )
p8  
ggsave(paste0(out_dir,'plot/sample_info/bank_sex_adnc_sample_information.pdf'), 
       plot = p8, 
       width = 8, 
       height = 10)

# protein missing rate
missing_rate_df <- data.frame(missing_rate = pro_raw$missing_rate)

p9 <- ggplot(missing_rate_df, aes(x = missing_rate)) +
  geom_histogram(binwidth = 0.05, color = "black", fill = "#2171B5", alpha = 0.7) +
  geom_text(stat = "bin", binwidth=0.05, aes(label = after_stat(count)), vjust = -0.5, size = 4) +
  labs(title = "Missing rate of protein",
       x = "Missing Rate",
       y = "Frequency") +
  theme_bw() + 
  theme(
    panel.grid = element_line(color = 'white'),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    plot.title = element_text(size = 20)
  )
p9
ggsave(paste0(out_dir,'plot/sample_info/protein_missing_rate.pdf'), 
       plot = p9, 
       width = 8, 
       height = 6)

# all sample boxplot
pro_50_norm <- as.data.frame(fread('processed_data/human_brain_pro_729_50_norm.txt'))
pro_box <- pro_50_norm[,colnames(pro_50_norm)%in%sample_adnc$id] %>%
  pivot_longer(
    cols = everything(),  # 选择所有列作为值列
    names_to = "Sample",   # 列名存入 "Sample"
    values_to = "Expression",  # 值存入 "Expression"
    values_drop_na = FALSE  # 是否删除 NA
  )
pro_box$adnc <- sample_adnc$adnc_num[match(pro_box$Sample, sample_adnc$id)]
pro_box$Disease <- sample_adnc$disease[match(pro_box$Sample, sample_adnc$id)]
pro_box <- pro_box[order(pro_box$Disease, pro_box$Sample), ]  # 按 Disease 先排序
pro_box$Sample <- factor(pro_box$Sample, levels = unique(pro_box$Sample))

hc_median_mean <- pro_box %>%
  filter(Disease == "HC") %>%
  group_by(Sample) %>%
  summarise(median_expr = median(Expression, na.rm = TRUE)) %>%
  summarise(mean_median_expr = mean(median_expr, na.rm = TRUE)) %>%
  pull(mean_median_expr)

p10 <- ggplot(pro_box, aes(x = reorder(Sample, Disease), y = Expression, fill = Disease)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +  # 箱线图，隐藏离群点
  geom_hline(yintercept = hc_median_mean, linetype = "dashed", color = "red", linewidth = 1) + # 添加水平线
  scale_fill_manual(values = c("HC" = "#C6DBEF", "ADNC-L" = "#6BAED6", "ADNC-M" = "#2171B5", "ADNC-H" = "#08306B")) + 
  coord_cartesian(ylim = c(quantile(pro_box$Expression, 0.01, na.rm = TRUE), 
                           quantile(pro_box$Expression, 0.99, na.rm = TRUE))) +
  labs(title = "",x = "",y = "") +
  theme_bw() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 25),
        panel.grid = element_line(color = 'white'))  # 倾斜 x 轴标签
p10
ggsave(paste0(out_dir,'plot/sample_info/all_sample_boxplot_norm.png'), 
       plot = p10, 
       width = 30, 
       height = 20,
       dpi = 300)

# disease type
disease_info <- sample_adnc[,c(1,36,52:58)]
other_disease <- pivot_longer(
  disease_info,
  cols = c(part, late, artag, lbd, cvd),
  names_to = "disease_type",
  values_to = "status"
)
other_disease <- other_disease[!is.na(other_disease$status),]
other_disease$Status <- ifelse(other_disease$status==0,'Control','Case')
p11 <- ggplot(other_disease, aes(x = disease_type, fill = Status)) +
  geom_bar(position = "dodge") +
  geom_text(stat = "count", 
            aes(label = after_stat(count)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 5) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(x = NULL, y = "Count", fill = "Status") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_line(color = 'white'),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(size = 20)
  )
p11
ggsave(paste0(out_dir,'plot/sample_info/other_disease.pdf'), 
       plot = p11, 
       width = 8, 
       height = 6)

# upset (ADNC CVD LATE Dementia et al)
ADNC_sample = sample_adnc$id[which(sample_adnc$adnc==1)]
dementia_sample = sample_adnc$id[which(sample_adnc$diag_dementia==1)]
cvd_sample = sample_adnc$id[which(sample_adnc$cvd==1)]
part_sample = sample_adnc$id[which(sample_adnc$part==1)]
late_sample = sample_adnc$id[which(sample_adnc$late==1)]
artag_sample = sample_adnc$id[which(sample_adnc$artag==1)]
sample_list <- list('ADNC'=ADNC_sample,'Dementia'=dementia_sample,
                    'CVD'=cvd_sample,'LATE'=late_sample,
                    'PART'=part_sample,'ARTAG'=artag_sample)

pdf(paste0(out_dir,'plot/sample_info/upset.pdf'), width = 10, height = 8)
upset(fromList(sample_list), 
      sets = c('ADNC','Dementia','CVD','LATE','PART','ARTAG'),
      point.size=2.5,  
      line.size=0.5, 
      mainbar.y.label="Count", 
      main.bar.color = "#2171B5",  
      sets.bar.color = "#C6DBEF",  
      sets.x.label="Numbers of sample",
      mb.ratio = c(0.7, 0.3), 
      order.by = "degree",
      decreasing = FALSE,
      text.scale = c(2, 1.8, 1.8, 1.8, 1.8, 1.8)
)
dev.off()
