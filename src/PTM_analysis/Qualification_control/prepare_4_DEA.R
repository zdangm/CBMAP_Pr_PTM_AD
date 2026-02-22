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
    library(ggrepel)
})
rm(list = ls())
conflicted::conflict_prefer_all("dplyr")
OUTLIER_LIST <- c("A2021CBB015")
if (!dir.exists(RES_SUB_DIR)) {
    dir.create(RES_SUB_DIR, recursive = T)
}

rep_code_2_id <-
    data.frame(
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


intensity_ms_idnttfd_info <- fread(
    paste0(RES_SUB_DIR, "/L0G0/MS_identified_information.txt")
)
intensity_ms_idnttfd_info.1 <- intensity_ms_idnttfd_info %>%
    select(starts_with("Intensity"))
dim(intensity_ms_idnttfd_info.1)
colnames(intensity_ms_idnttfd_info.1) <- gsub(
    "Intensity ",
    "",
    colnames(intensity_ms_idnttfd_info.1)
)
intensity_ms_idnttfd_info.1 <- intensity_ms_idnttfd_info.1 %>%
    select(-starts_with("MIX"))
dim(intensity_ms_idnttfd_info.1)
intensity_ms_idnttfd_info <- cbind(
    intensity_ms_idnttfd_info[, 1:7],
    intensity_ms_idnttfd_info.1
)
dim(intensity_ms_idnttfd_info)
if (F) {
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
    rep_code_2_id
    tmp <- intensity_ms_idnttfd_info %>%
        select(any_of(rep_code_2_id$repeat_code), any_of(rep_code_2_id$id))
    tmp.cv <- apply(rep_code_2_id, 1, function(x) {
        if ((x[1] %in% colnames(tmp)) & (x[2] %in% colnames(tmp))) {} else {
            return(NA)
        }
        tmp.x <- tmp %>%
            select(all_of(x[1]), all_of(x[2])) %>%
            as.matrix() %>%
            as.vector()
        tmp.x <- log2(tmp.x)
        mean.x <- mean(tmp.x, na.rm = TRUE)
        sd.x <- sd(tmp.x, na.rm = TRUE)
        cv.x <- sd.x / mean.x
        cv.x
    })
    tmp.cv
    tmp.cv.df <- cbind(rep_code_2_id, cv = tmp.cv)
    tmp.cv.df
    write.xlsx(
        tmp.cv.df,
        file = paste0(RES_SUB_DIR, "repeat_code_cv.xlsx"),
        rowNames = FALSE
    )
}


rm(intensity_ms_idnttfd_info.1)
gc()

ms_idnttfd_info <- fread(
    paste0(RES_SUB_DIR, "/L0G1/MS_identified_information.txt")
)
ms_idnttfd_info <- ms_idnttfd_info %>% select(-starts_with("Intensity"))
ms_idnttfd_info <- ms_idnttfd_info %>% select(-starts_with("MIX"))

ms_info.ls <- list(
    intensity = intensity_ms_idnttfd_info,
    no_intensity = ms_idnttfd_info
)
rm(ms_idnttfd_info, intensity_ms_idnttfd_info)
gc()


ms_info.ls <- lapply(ms_info.ls, function(x) {
    x <- x %>% rename(PTB087 = repeat2)
    x
})
lapply(ms_info.ls, dim)

ms_info.ls <- lapply(ms_info.ls, function(x) {
    # x$identifier <- paste(x$`Protein accession`, x$Position, x$`Amino acid`, sep = "_")
    x <- x %>%
        group_by(`Gene name`, Position, `Amino acid`) %>%
        filter(!str_detect(`Protein accession`, "_")) %>%
        ungroup()
    x$identifier <- paste(
        x$`Protein accession`,
        x$Position,
        x$`Amino acid`,
        sep = "_"
    )
    x$identifier <- paste(x$`Gene name`, x$identifier, sep = ".")
    x <- x %>% select(identifier, everything())
    x
})


sample_list_765 <- read.xlsx(
    "sample_meta_info/pr_measured_ZJU_765_id.xlsx",
    sheet = 1
)
sample_list_765 <- rbind(sample_list_765, c("PTB087", "PTB087", "PTB087"))


ms_info.723.ls <- lapply(ms_info.ls, function(x) {
    x.730 <- x %>%
        select(all_of(
            sample_list_765[
                !is.na(sample_list_765$CBMAP_profile_1187),
            ]$jingjie_ID
        ))
    x.723 <- x.730 %>% select(-all_of(OUTLIER_LIST))
    x.723 <- cbind(x[, 1:8], x.723)
    x.723
})


rm(ms_info.ls)
gc()
missing_rate.ls <- lapply(ms_info.723.ls, function(x) {
    missing_rate <- rowMeans(is.na(x[, -c(1:8)]))
    missing_rate
})
missing_df <- data.frame(missing_rate = missing_rate.ls[[1]])

p.Phos.miss.rate <- ggplot(missing_df, aes(x = missing_rate)) +
    geom_histogram(
        binwidth = 0.1,
        boundary = 0,
        fill = "#0680be",
        closed = "left",
        color = "black"
    ) + # 设置柱宽
    geom_text(
        stat = "bin",
        aes(label = after_stat(count)),
        vjust = -0.5,
        binwidth = 0.1,
        size = 4,
        boundary = 0
    ) +
    labs(
        title = "Distribution of Phos-site Missing Rates",
        x = "Missing Rate",
        y = "Frequency"
    ) +
    theme_gray() +
    theme(
        axis.text.x = element_text(size = 16),
        panel.grid.minor = element_blank(),
        panel.grid = element_line(color = "white"),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.8
        ),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = 16)
    ) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1))

pdf(file = "missing_rate.pdf", width = 8, height = 6)
print(p.Phos.miss.rate)
dev.off()


qs_savem(
    ms_info.723.ls,
    file = paste0(
        RES_SUB_DIR,
        "intensity_and_proteomic_scaled_ms_idnttfd_info.729.QS2"
    )
)
lapply(ms_info.723.ls, dim)
ms_info.723.missing.let0.5.ls <- mapply(
    function(x, y) {
        x[y < 0.5, ]
    },
    ms_info.723.ls,
    missing_rate.ls,
    SIMPLIFY = F
)

rm(ms_info.723.ls)
gc()


rm(missing_rate.ls)
gc()
qs_savem(
    ms_info.723.missing.let0.5.ls,
    file = paste0(
        RES_SUB_DIR,
        "intensity_and_proteomic_scaled_ms_idnttfd_info.723.missing.let0.5.QS2"
    )
)


#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------#
log1p_raw.ls <- lapply(ms_info.723.missing.let0.5.ls, function(x) {
    log1p_raw <- log2(x[, -c(1:8)] + 1)
    log1p_raw <- as.matrix(log1p_raw)
    rownames(log1p_raw) <- x$identifier
    log1p_raw <- t(log1p_raw)
    log1p_raw <- data.frame(log1p_raw, check.names = F)
    log1p_raw
})

rm(ms_info.723.missing.let0.5.ls)
gc()


#---------------------------------#
#-------prepare covariates--------#
#---------------------------------#
batch_info <- read.xlsx(
    "sample_meta_info/【样本明细】966标+21标MIX（共987标）人脑组织4D fastDIA蛋白组、磷酸化、乙酰化、泛素化项目正式样本明细20241121.xlsx",
    sheet = 3
)

batch_info[match("repeat2", batch_info$`分析名`), "分析名"] <- "PTB087"
batch_info <- batch_info[
    match(rownames(log1p_raw.ls[[1]]), batch_info$分析名),
    c("分析名", "样本标签", "IP批次", "酶解批次")
]


sample_info <- fread(
    "./sample_meta_info/sample_info_maintained_by_ZLab/CBMAP_sample_info_1187_final_20250523.csv"
)


sample_info$id <- if_else(
    sample_info$bank == "zju",
    paste0("A", sample_info$id),
    sample_info$id
)
sample_info <- sample_info[match(rownames(log1p_raw.ls[[1]]), sample_info$id), ]

covariates <- batch_info %>% left_join(sample_info, by = c("分析名" = "id"))
colnames(covariates)
colnames(covariates)[1:4] <- c("id", "sample_label", "IP_batch", "digest_batch")
colnames(covariates) <- gsub("-", "_", colnames(covariates))
colnames(covariates) <- gsub(" ", "_", colnames(covariates))
colnames(covariates) <- gsub("\\.", "_", colnames(covariates))
colnames(covariates)
dim(covariates)
rm(sample_info, batch_info)
gc()


covariates <- covariates[match(rownames(log1p_raw.ls[[1]]), covariates$id), ]
log1p_raw.all.info.ls <- lapply(log1p_raw.ls, function(x) {
    cbind(x, covariates)
})


#---------------------------------------------------------------------------------------------------------------#
log1p_raw.all.info.ls <- lapply(log1p_raw.all.info.ls, function(x) {
    x <- x %>%
        mutate(
            ADNC_LMH = case_when(
                is.na(ADNC) ~ NA,
                ADNC == 0 ~ 0,
                ADNC_L == 1 ~ 1,
                ADNC_M == 1 ~ 2,
                ADNC_H == 1 ~ 3
            )
        )
    x <- x %>%
        mutate(
            hc = case_when(
                is.na(ADNC) ~ NA,
                ADNC == 0 ~ 1,
                ADNC == 1 ~ 0
            )
        )
    x$sex_male <- factor(x$sex_male, ordered = FALSE)
    x
})
lapply(log1p_raw.all.info.ls, dim)

rm(log1p_raw.ls)
gc()
qs_savem(
    log1p_raw.all.info.ls,
    file = paste0(RES_SUB_DIR, "log2_raw.729.all.info.ls.QS2")
)


qs_readm(paste0(RES_SUB_DIR, "log2_raw.729.all.info.ls.QS2"))
if (PTMs[i] == "ace") {
    NUM_COV <- 56
    log1p_raw.all.info.ls <- lapply(log1p_raw.all.info.ls, function(x) {
        x.probes <- x[, seq.int(1, ncol(x) - NUM_COV)]
        x.covs <- x[, (ncol(x) - NUM_COV + 1):ncol(x)]
        x.probes <- 2^x.probes - 1
        x <- cbind(x.probes, x.covs)
        x
    })

    ppline <- "log_impute"

    pipeline_steps <- strsplit(ppline, "_")[[1]]

    for (step in pipeline_steps) {
        func <- preprocess_in_DEA[[step]]
        if (is.null(func)) {
            stop(paste("wrong: ", step))
        }
        log1p_raw.all.info.ls <- lapply(log1p_raw.all.info.ls, func, NUM_COV)
    }
    log1p_raw.all.info.ls.bak <- log1p_raw.all.info.ls
    log1p_raw.all.info.ls <- log1p_raw.all.info.ls.bak

    library(sva)
    log1p_raw.all.info.ls <- lapply(log1p_raw.all.info.ls, function(x) {
        x <- filter(x, !is.na(ADNC_LMH))
        x.mat <- x[, seq.int(1, ncol(x) - NUM_COV)]
        x.cov <- x[, -seq.int(1, ncol(x) - NUM_COV)]
        ip_batch <- if_else(x.cov$IP_batch == "B18", 1, 0)
        ip_batch <- factor(ip_batch)
        x.cov$ADNC_LMH <- factor(x.cov$ADNC_LMH, ordered = T)

        for_mod <- model.matrix(
            as.formula("~ ADNC_LMH + age + sex_male"),
            x.cov
        )
        x.cov$ADNC_LMH <- as.numeric(x.cov$ADNC_LMH) - 1
        x.mat.t <- t(x.mat)
        x.adj.t <- ComBat(dat = x.mat.t, batch = ip_batch, mod = for_mod)
        x.adj <- t(x.adj.t)
        x <- cbind(x.adj, x.cov)
        x
    })
}

ppline <- "scale_INT"

pipeline_steps <- strsplit(ppline, "_")[[1]]
for (step in pipeline_steps) {
    func <- preprocess_in_DEA[[step]]
    if (is.null(func)) {
        stop(paste("❌ 无效的处理步骤：", step))
    }
    log1p_raw.all.info.ls <- lapply(log1p_raw.all.info.ls, func, NUM_COV)
}
lapply(log1p_raw.all.info.ls, dim)
if (PTM[i] == "ace") {
    qs_savem(
        log1p_raw.all.info.ls,
        file = paste0(
            RES_SUB_DIR,
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
        )
    )
} else {
    qs_savem(
        log1p_raw.all.info.ls,
        file = paste0(
            RES_SUB_DIR,
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.QS2"
        )
    )
}
