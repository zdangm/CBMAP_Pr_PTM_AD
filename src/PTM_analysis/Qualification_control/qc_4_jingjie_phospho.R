suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(openxlsx)
  library(magrittr)
  library(patchwork)
  library(limma)
  library(UpSetR)
  library(gridExtra)
  library(qs2)
  library(optparse)
  library(purrr)
})
rm(list = ls())
conflicted::conflict_prefer_all("dplyr")

OUTLIER_LIST <- c("A2021CBB015")
rep_code_2_id <- data.frame(
  repeat_code = paste0("repeat", 1:7),
  id = c(
    "PTB425",
    "PTB087",
    "XYA20221024",
    "XYA20191227",
    "A2024CBB001",
    "A2023CBB072",
    "A2022CBB082"
  )
)
PTMs <- c("phos", "ubiq", "ace")

option_list <- list(
  make_option(
    c("--params"),
    type = "integer",
    default = 1,
    help = "ptm",
    metavar = "number"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
i <- opt$params
RES_SUB_DIR <- sprintf("../data/%s/", PTMs[i])

ms_idnttfd_info <- fread(
  paste0(
    RES_SUB_DIR,
    "/data/XB07045B4DPST_0311/XB07045B4DPST_mix+sample/L0G0/MS_identified_information.txt"
  )
)
dim(ms_idnttfd_info)
intensity_ms_idnttfd_info <- ms_idnttfd_info %>%
  select(starts_with("Intensity"))
dim(intensity_ms_idnttfd_info)
# 111213    987
colnames(intensity_ms_idnttfd_info) <- gsub(
  "Intensity ",
  "",
  colnames(intensity_ms_idnttfd_info)
)
intensity_ms_idnttfd_info <- cbind(
  ms_idnttfd_info[, 1:7],
  intensity_ms_idnttfd_info
)
intensity_ms_idnttfd_info$identifier <- paste(
  intensity_ms_idnttfd_info$`Protein accession`,
  intensity_ms_idnttfd_info$Position,
  intensity_ms_idnttfd_info$`Amino acid`,
  sep = "_"
)
intensity_ms_idnttfd_info <- intensity_ms_idnttfd_info %>%
  select(identifier, starts_with("MIX"), everything())
dim(intensity_ms_idnttfd_info)


scld_intensity_ms_idnttfd_info <- ms_idnttfd_info %>%
  select(!starts_with("Intensity"))
scld_intensity_ms_idnttfd_info %>% dim()
scld_intensity_ms_idnttfd_info <- scld_intensity_ms_idnttfd_info[, -c(1:7)]

all.equal(
  is.na(intensity_ms_idnttfd_info),
  is.na(scld_intensity_ms_idnttfd_info)
)
rm(ms_idnttfd_info, scld_intensity_ms_idnttfd_info)
# F
gc()


sample_list_765 <- read.xlsx(
  "sample_meta_info/pr_measured_ZJU_765_id.xlsx",
  sheet = 1
)
sample_list_765 %>% dim()

sample_infor <- read.csv(
  "sample_information/final/CBMAP_sample_info_1187_final_20250523.csv"
)
sample_infor <- sample_infor %>% select(-disease, -diagnosis)
sample_infor %>% dim()
sample_infor$id <- if_else(
  sample_infor$bank == "zju",
  paste0("A", sample_infor$id),
  sample_infor$id
)
sample_infor <- sample_infor[sample_infor$id %in% sample_list_765$jingjie_ID, ]
dim(sample_infor)


batch_infor <- read.xlsx(
  "sample_meta_info/【样本明细】966标+21标MIX（共987标）人脑组织4D fastDIA蛋白组、磷酸化、乙酰化、泛素化项目正式样本明细20241121.xlsx",
  sheet = 3
)
batch_infor <- batch_infor %>% select(`分析名`, `酶解批次`, `IP批次`)
batch_infor <- batch_infor %>%
  left_join(rep_code_2_id, by = c("分析名" = "repeat_code"))
batch_infor$id <- if_else(
  is.na(batch_infor$id),
  batch_infor$分析名,
  batch_infor$id
)
dim(batch_infor)
which(duplicated(batch_infor$id))
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------#
#-------------------------------------------Excluding outliers---------------------------------------------------#
#-----------------------------using shared ub-sites across all samples to perform PCA----------------------------#
#----------------------------------------------------------------------------------------------------------------#
intensity_ms_idnttfd_info.1 <- intensity_ms_idnttfd_info %>%
  dplyr::select(
    all_of(sample_list_765$jingjie_ID),
    starts_with("repeat"),
    starts_with("MIX")
  )
dim(intensity_ms_idnttfd_info.1)
intensity_ms_idnttfd_info.1 <- cbind(
  intensity_ms_idnttfd_info[, 1],
  intensity_ms_idnttfd_info.1
)
row_miss_rate <- rowMeans(is.na(intensity_ms_idnttfd_info.1[, -1]))


intensity_ms_idnttfd_info.2 <- intensity_ms_idnttfd_info.1

intensity_ms_idnttfd_info.2.no.missing <- intensity_ms_idnttfd_info.2[
  row_miss_rate == 0,
]
dim(intensity_ms_idnttfd_info.2.no.missing)
intensity_ms_idnttfd_info.2.no.missing <- as.matrix(
  intensity_ms_idnttfd_info.2.no.missing
)
rownames(
  intensity_ms_idnttfd_info.2.no.missing
) <- intensity_ms_idnttfd_info.2.no.missing[, 1]
intensity_ms_idnttfd_info.2.no.missing <- intensity_ms_idnttfd_info.2.no.missing[,
  -1
]
intensity_ms_idnttfd_info.2.no.missing.t <- t(
  intensity_ms_idnttfd_info.2.no.missing
)
intensity_ms_idnttfd_info.2.no.missing.t <- apply(
  intensity_ms_idnttfd_info.2.no.missing.t,
  2,
  as.numeric
)
pca <- prcomp(intensity_ms_idnttfd_info.2.no.missing.t, scale. = T)


pca$x %>% dim()
pca.x <- as.data.frame(pca$x)
pca.x$id <- colnames(intensity_ms_idnttfd_info.2.no.missing)
pca.x <- pca.x %>% select(id, everything())
dim(pca.x)

rep_code_2_id
rep_code_2_id$code <- str_extract(rep_code_2_id$repeat_code, "\\d")
rep_code_2_id$code <- paste0("original", rep_code_2_id$code)
rep_code_2_id <- rep_code_2_id %>%
  filter(repeat_code %in% pca.x$id, id %in% pca.x$id)


pca.x <- pca.x %>% left_join(rep_code_2_id[, -1], by = "id")
table(pca.x$code, useNA = "always")
pca.x$code <- if_else(
  is.na(pca.x$code),
  if_else(
    str_detect(pca.x$id, "repeat") & pca.x$id %in% rep_code_2_id$repeat_code,
    pca.x$id,
    NA
  ),
  pca.x$code
)
table(pca.x$code, useNA = "always")
pca.x %>% filter(!is.na(code)) %>% select(id, code)
rep_code_2_id
pca.x <- pca.x %>%
  mutate(
    state = case_when(
      is.na(code) ~ "other",
      str_detect(code, "^repeat") ~ "repeat",
      str_detect(code, "^original") ~ "original"
    )
  )
pca.x$code <- if_else(
  is.na(pca.x$code),
  "other",
  str_extract(pca.x$code, "\\d")
)
pca.x$code <- pca.x$code %>% dense_rank() %>% as.character()
pca.x$code <- if_else(pca.x$code == "6", "other", paste0("pair", pca.x$code))
pca.x <- pca.x %>% select(id, code, state, everything())
pca.x$state <- factor(
  pca.x$state,
  levels = c("other", "repeat", "original"),
  ordered = T
)


pca.repsamples <- ggplot() +
  geom_point(
    data = subset(pca.x, state == "other"),
    aes(x = PC1, y = PC2, shape = state, color = code),
    size = 2,
    alpha = 1
  ) +
  geom_point(
    data = subset(pca.x, state != "other"),
    aes(x = PC1, y = PC2, shape = state, color = code),
    size = 2,
    alpha = 1
  ) +
  scale_color_manual(
    values = c(
      "pair1" = "#A6CEE3",
      "pair2" = "#B2DF8A",
      "pair3" = "#FDBF6F",
      "pair4" = "#CAB2D6",
      "pair5" = "#FF9999",
      "other" = "lightgray"
    )
  ) +
  stat_ellipse(
    data = pca.x,
    aes(x = PC1, y = PC2, color = code),
    level = 0.95
  ) +
  theme_bw() +
  theme(
    panel.grid = element_line(color = "white"),
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 14)
  )
pca.repsamples


pca.x$Type <- str_detect(pca.x$id, "MIX")
pca.x$Type <- if_else(pca.x$Type, "MIX", "Non-MIX")
pca.x <- pca.x %>% select(id, Type, everything()) %>% arrange(Type)

pca.mix <- ggplot() +
  geom_point(
    data = subset(pca.x, Type == "Non-MIX"),
    aes(x = PC1, y = PC2, color = Type),
    size = 2,
    alpha = 0.9
  ) +
  geom_point(
    data = subset(pca.x, Type == "MIX"),
    aes(x = PC1, y = PC2, color = Type),
    size = 2,
    alpha = 0.9
  ) +
  scale_color_manual(values = c("MIX" = "#FF9999", "Non-MIX" = "lightgray")) +
  stat_ellipse(
    data = pca.x,
    aes(x = PC1, y = PC2, color = Type),
    level = 0.95
  ) +
  theme_bw() +
  # ggtitle("MIX Sample") +
  theme(
    panel.grid = element_line(color = "white"),
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
pca.mix


pca.x.sex <- pca.x %>%
  inner_join(sample_infor, by = c("id" = "id")) %>%
  select(PC1, PC2, sex_male)
pca.x.sex %>% dim()
pca.x.sex$sex_male <- factor(pca.x.sex$sex_male, labels = c("female", "male"))
pca.sex <- ggplot(pca.x.sex, aes(PC1, PC2, color = sex_male)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  # ggtitle("sex") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 14),
    # legend.position = "none",
    panel.grid = element_line(color = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 14)
  )

pca.sex

pca.x.bank <- pca.x %>%
  inner_join(sample_infor, by = c("id" = "id")) %>%
  select(PC1, PC2, bank)
pca.x.bank$Bank <- as.factor(pca.x.bank$bank)
dim(pca.x.bank)
pca.bank <- ggplot(pca.x.bank, aes(PC1, PC2, color = Bank)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  # ggtitle("Bank") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 14),
    # legend.position = "none",
    panel.grid = element_line(color = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 14)
  )
pca.bank


pca.x.enzyme <- pca.x %>%
  inner_join(batch_infor, by = c("id" = "分析名")) %>%
  select(PC1, PC2, `酶解批次`)
dim(pca.x.enzyme)
pca.x.enzyme$enzyme_batch <- as.factor(pca.x.enzyme$酶解批次)
pca.enzyme <- ggplot(pca.x.enzyme, aes(PC1, PC2, color = enzyme_batch)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  # ggtitle("Enzyme batch") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 14),
    # legend.position = "none",
    panel.grid = element_line(color = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 14)
  )
pca.enzyme


pca.x.IP <- pca.x %>%
  inner_join(batch_infor, by = c("id" = "分析名")) %>%
  select(PC1, PC2, `IP批次`)
pca.x.IP$ip_batch <- factor(
  pca.x.IP$IP批次,
  levels = gtools::mixedsort(unique(pca.x.IP$IP批次))
)
dim(pca.x.IP)
pca.IP <- ggplot(pca.x.IP, aes(PC1, PC2, color = ip_batch)) +
  geom_point() +
  stat_ellipse(level = 0.95) +
  # ggtitle("IP batch") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_text(size = 14),
    # legend.position = "none",
    panel.grid = element_line(color = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = 14)
  ) +
  guides(color = guide_legend(ncol = 2))
pca.IP


library(patchwork)
pdf("post_mut_pca.pdf", width = 10, height = 10)
pca.repsamples +
  pca.mix +
  pca.sex +
  pca.bank +
  pca.enzyme +
  pca.IP +
  plot_layout(ncol = 2, nrow = 3)
dev.off()
p.ls <- list(pca.repsamples, pca.mix, pca.sex, pca.bank, pca.enzyme, pca.IP)
pdf("pca.pdf", width = 7, height = 5)
for (p in seq_along(p.ls)) {
  print(p.ls[[p]])
}
dev.off()
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------#
#-------------------------------------------check consistency----------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
dim(intensity_ms_idnttfd_info.1)
intensity_ms_idnttfd_info.1 <- intensity_ms_idnttfd_info.1 %>%
  select(-starts_with("MIX"))


compute_corr_matrix <- function(df, replicate_map) {
  sample_cols <- colnames(df)[-1]

  replicate_samples <- replicate_map[, 1]
  original_samples <- replicate_map[, 2]

  other_samples <- setdiff(sample_cols, replicate_samples)
  other_samples <- setdiff(other_samples, original_samples)

  result_matrix <- matrix(
    NA,
    nrow = length(replicate_samples),
    ncol = length(original_samples) + length(other_samples),
    dimnames = list(replicate_samples, c(original_samples, other_samples))
  )

  for (i in seq_along(replicate_samples)) {
    replicate_sample <- replicate_samples[i]
    original_sample <- original_samples[i]

    replicate_data <- df[[replicate_sample]]
    original_data <- df[[original_sample]]

    corr <- compute_correlation(replicate_data, original_data)
    result_matrix[i, original_sample] <- corr

    for (other_sample in other_samples) {
      other_data <- df[[other_sample]]
      corr_other <- compute_correlation(replicate_data, other_data)
      result_matrix[i, other_sample] <- corr_other
    }
  }

  return(result_matrix)
}

compute_correlation <- function(x, y) {
  valid_idx <- !is.na(x) & !is.na(y)
  x_valid <- x[valid_idx]
  y_valid <- y[valid_idx]

  if (length(x_valid) < 2 || length(y_valid) < 2) {
    return(NA)
  }

  lm_model <- lm(y_valid ~ x_valid)
  p_value <- summary(lm_model)$coefficients[2, 4]

  if (p_value < 0.05) {
    corr <- cor(x_valid, y_valid, method = "pearson", use = "complete.obs")
  } else {
    corr <- cor(x_valid, y_valid, method = "spearman", use = "complete.obs")
  }

  return(corr)
}


plot_boxplot <- function(
  matrix_data,
  main = "boxplot of corr",
  range = c(0, 1)
) {
  n <- nrow(matrix_data)

  diag_elements <- diag(as.matrix(matrix_data[, 1:n], nrow = n))

  off_diag_elements <- as.vector(as.matrix(matrix_data[,
    (n + 1):ncol(matrix_data)
  ]))

  diag_elements <- diag_elements[!is.na(diag_elements)]
  off_diag_elements <- off_diag_elements[!is.na(off_diag_elements)]

  data <- data.frame(
    Group = rep(
      c("repeat-original", "repeat-others"),
      times = c(length(diag_elements), length(off_diag_elements))
    ),
    Value = c(diag_elements, off_diag_elements)
  )

  p <- ggplot(data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot() +
    # geom_jitter(aes(color = Group), alpha = 0.3, width = 0.2, size = 1) +
    theme_minimal() +
    labs(
      title = "",
      x = "Group",
      y = "Correlation Value"
    ) +
    coord_cartesian(ylim = range) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.text.x = element_blank(),
      panel.grid = element_line(color = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title = element_blank(),
      legend.text = element_text(size = 20),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.title = element_blank()
    )
  p
}


df1c <- intensity_ms_idnttfd_info.1[, 1]
logdf <- log2(1 + intensity_ms_idnttfd_info.1[, -1])
logdf <- cbind(df1c, logdf)
log.corr.mat <- compute_corr_matrix(logdf, rep_code_2_id)
log.corr.mat[, 1:10]
log.corr.boxplot <- plot_boxplot(log.corr.mat, range = c(0.7, 1))
pdf("post_mut_log_rep_corr.pdf", width = 8, height = 5)
log.corr.boxplot
dev.off()
