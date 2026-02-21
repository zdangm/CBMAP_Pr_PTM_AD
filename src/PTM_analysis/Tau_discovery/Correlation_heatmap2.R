suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(openxlsx)
    library(magrittr)
    library(patchwork)
    library(limma)
    library(qs2)
    library(purrr)
    library(ggrepel)
    library(impute)
    library(clusterProfiler)
    library(circlize)
    library(ComplexHeatmap)
    library(dplyr)
})
rm(list = ls())
NUM_COV <- 56
conflicted::conflict_prefer_all("dplyr")
conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(base::setdiff)
source(
    "helper_func.R"
)

RES_SUB_DIR.ls <- list(
    phospho = "phos/",
    ubiq = "ubiq/",
    ace = "ace/"
)


is_tau <- T
#--------------------------------------------------Begin---------------------------------------------------------#
MAPT.ls.1 <- list(intensity = NULL, no_intensity = NULL)
for (idx in seq_along(RES_SUB_DIR.ls)) {
    RES_SUB_DIR <- RES_SUB_DIR.ls[[idx]]
    ptm <- names(RES_SUB_DIR.ls)[idx]
    if (ptm == "ace") {
        qs_readm(paste0(
            RES_SUB_DIR,
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
        ))
    } else {
        qs_readm(paste0(
            RES_SUB_DIR,
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.QS2"
        ))
    }
    # qs_readm(paste0(RES_SUB_DIR, "log2_raw.729.all.info.ls.QS2"))

    for (j in seq_along(log1p_raw.all.info.ls)) {
        x <- log1p_raw.all.info.ls[[j]]
        # x <- myimputeknn(x, NUM_COV = 56)
        # x <- mymedian(x, NUM_COV = 56)
        # x <- myINT(x, NUM_COV = 56)
        if (is_tau) {
            x <- x %>% select(matches("^MAPT\\.")) # for Tau
        } else {
            x <- x %>% select(matches("^APP\\.")) # for APP
        }
        colnames(x) <- paste0(ptm, ".", colnames(x))
        log1p_raw.all.info.ls[[j]] <- x
    }

    for (j in seq_along(MAPT.ls.1)) {
        df <- MAPT.ls.1[[j]]
        if (is.null(df)) {
            df <- log1p_raw.all.info.ls[[j]]
        } else {
            df <- df[
                match(rownames(log1p_raw.all.info.ls[[j]]), rownames(df)),
            ]
            df <- cbind(df, log1p_raw.all.info.ls[[j]])
        }
        MAPT.ls.1[[j]] <- df
    }
}


#---------------------------------#
#-------prepare covariates--------#
#---------------------------------#
batch_info <- read.xlsx(
    "/share/home/lik/Scripts/likun/sample_meta_info/【样本明细】966标+21标MIX（共987标）人脑组织4D fastDIA蛋白组、磷酸化、乙酰化、泛素化项目正式样本明细20241121.xlsx",
    sheet = 3
)
dim(batch_info)
# 987 13
"PTB087" %in% batch_info$`分析名`
# FALSE
"repeat2" %in% batch_info$`分析名`
# TRUE
batch_info[match("repeat2", batch_info$`分析名`), "分析名"] <- "PTB087"
batch_info <- batch_info[
    match(rownames(MAPT.ls.1[[1]]), batch_info$分析名),
    c("分析名", "样本标签", "IP批次", "酶解批次")
]

sample_info <- fread(
    "/data/shared_data/China_Brain_MultiOmics/sample_information/final/CBMAP_sample_info_1187_final_20250523.csv"
)
sample_info$id <- if_else(
    sample_info$bank == "zju",
    paste0("A", sample_info$id),
    sample_info$id
)
sample_info <- sample_info[match(rownames(MAPT.ls.1[[1]]), sample_info$id), ]

covariates <- batch_info %>% left_join(sample_info, by = c("分析名" = "id"))
colnames(covariates)
colnames(covariates)[1:4] <- c("id", "sample_label", "IP_batch", "digest_batch")
colnames(covariates) <- gsub("-", "_", colnames(covariates))
colnames(covariates) <- gsub(" ", "_", colnames(covariates))
colnames(covariates) <- gsub("\\.", "_", colnames(covariates))
colnames(covariates)

rm(sample_info, batch_info)
gc()


#----------------------------------------------------------#
#-------Using limma and glm to perform DEA-----------------#
#----------------------------------------------------------#

covariates <- covariates[match(rownames(MAPT.ls.1[[1]]), covariates$id), ]
identical(rownames(MAPT.ls.1[[1]]), covariates$id)
identical(rownames(MAPT.ls.1[[1]]), rownames(MAPT.ls.1[[2]]))
# TRUE
MAPT.ls.2 <- lapply(MAPT.ls.1, function(x) {
    cbind(x, covariates)
})

df.hc.ls <- lapply(MAPT.ls.2, function(x) {
    df.hc <- x %>%
        filter(
            ADNC == 0,
            LBD == 0,
            ARTAG_Gray_matter == 0,
            ARTAG_Perivascular == 0,
            ARTAG_Subpial == 0,
            OTHER == 0,
            diag_dementia == 0,
            diag_PD == 0,
            diag_ALS == 0,
            diag_SCZ == 0,
            diag_epilepsy == 0,
            diag_others == 0
        )
    df.hc
})
lapply(df.hc.ls, dim)
# $intensity
# [1]   198 181

# $no_intensity
# [1]   198 181
df.ADNC.case.ls <- lapply(MAPT.ls.2, function(x) {
    df.ADNC.case <- x %>%
        filter(
            ADNC == 1,
            OTHER == 0,
            diag_ALS == 0,
            diag_SCZ == 0,
            diag_epilepsy == 0,
            diag_others == 0
        )
    df.ADNC.case
})
lapply(df.ADNC.case.ls, dim)
# $intensity
# [1]   252 181

# $no_intensity
# [1]   252 181
df.ADNC.ls <- mapply(
    function(x, y) {
        rbind(x, y)
    },
    df.hc.ls,
    df.ADNC.case.ls,
    SIMPLIFY = F
)
lapply(df.ADNC.ls, dim)
# $intensity
# [1] 450 181

# $no_intensity
# [1] 450 181

df.ADNC.ls <- lapply(df.ADNC.ls, function(x) {
    x <- x %>%
        mutate(
            ADNC_LMH = case_when(
                is.na(ADNC) ~ NA,
                ADNC == 0 ~ "HC",
                ADNC_L == 1 ~ "Low",
                ADNC_M == 1 ~ "Moderate",
                ADNC_H == 1 ~ "High"
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


lapply(df.ADNC.ls, dim)
sig.sites.ls <- list(intensity = NULL, no_intensity = NULL)
for (i in seq_along(sig.sites.ls)) {
    for (j in seq_along(RES_SUB_DIR.ls)) {
        all_sites.dir <- dirname(RES_SUB_DIR.ls[[j]])
        ptm <- names(RES_SUB_DIR.ls)[j]
        if (ptm == "ace") {
            all_sites.dir <- file.path(
                all_sites.dir,
                "version_after_2025_06_30/log_impute_scale_INT/All.sites.limma.QS2"
            )
        } else {
            all_sites.dir <- file.path(
                all_sites.dir,
                "version_after_2025_05_23/log_impute_scale_INT/All.sites.limma.QS2"
            )
        }
        qs_readm(all_sites.dir)
        All.sites <- lapply(All.sites, function(x) {
            x <- x %>%
                filter(
                    ADNC.LMH.num.P.Value < 0.05 | braak.num.P.Value < 0.05
                ) %>%
                select(
                    ADNC.LMH.num.P.Value,
                    braak.num.P.Value
                ) %>%
                rownames_to_column("identifier") %>%
                filter(str_detect(identifier, "^MAPT\\."))
            x$identifier <- paste(
                names(RES_SUB_DIR.ls)[j],
                x$identifier,
                sep = "."
            )
            x
        })
        All.sites <- All.sites[[i]]
        sig.sites.ls[[i]] <- rbind(
            sig.sites.ls[[i]],
            All.sites[, "identifier", drop = FALSE]
        )
    }
}
lapply(sig.sites.ls, dim)
lapply(df.ADNC.ls, dim)
colnames.ls <- list(intensity = NULL, no_intensity = NULL)
for (i in seq_along(df.ADNC.ls)) {
    df <- df.ADNC.ls[[i]]
    df <- df[, seq.int(1, ncol(df) - NUM_COV)]
    colnames.ls[[i]] <- data.frame(
        identifier = colnames(df),
        loc_on_MAPT = as.numeric(str_split_i(colnames(df), "_", 2)),
        aa = str_split_i(colnames(df), "_", 3),
        ptm = str_split_i(colnames(df), "\\.", 1)
    )
    colnames.ls[[i]]$loc_on_MAPT <- sapply(
        colnames.ls[[i]]$loc_on_MAPT,
        from_longest_isof_to_clinical_isof
    )
    colnames.ls[[i]] <- colnames.ls[[i]] %>%
        filter(!is.na(loc_on_MAPT), loc_on_MAPT != "not found on MAPT-8")
    colnames.ls[[i]]$loc_on_MAPT <- str_split_i(
        colnames.ls[[i]]$loc_on_MAPT,
        "_",
        2
    ) %>%
        as.numeric()
    colnames.ls[[i]] <- colnames.ls[[i]] %>% arrange(loc_on_MAPT)
    colnames.ls[[i]]$label <- paste0(
        colnames.ls[[i]]$ptm,
        ".",
        colnames.ls[[i]]$aa,
        colnames.ls[[i]]$loc_on_MAPT
    )
    colnames.ls[[i]]$is_sig <- if_else(
        colnames.ls[[i]]$identifier %in% sig.sites.ls[[i]]$identifier,
        TRUE,
        FALSE
    )
}


htmap.ls <- list(intensity = NULL, no_intensity = NULL)
for (i in seq_along(df.ADNC.ls)) {
    df <- df.ADNC.ls[[i]]
    df <- df[, seq.int(1, ncol(df) - NUM_COV)]
    # df <- df[, match(colnames.ls[[i]][colnames.ls[[i]]$is_sig == TRUE, ]$identifier, colnames(df))]
    df <- df[, match(colnames.ls[[i]]$identifier, colnames(df))]
    colnames.ls[[i]] <- colnames.ls[[i]] %>%
        mutate(
            domain = case_when(
                # loc_on_MAPT >= 1 & loc_on_MAPT <= 44 ~ "N-terminal",
                # loc_on_MAPT >= 45 & loc_on_MAPT <= 73 ~ "1N",
                # loc_on_MAPT >= 74 & loc_on_MAPT <= 102 ~ "2N",
                # loc_on_MAPT >= 103 & loc_on_MAPT <= 150 ~ "Linker (2N-PRR)",
                # loc_on_MAPT >= 151 & loc_on_MAPT <= 197 ~ "Proline-rich-1",
                # loc_on_MAPT >= 198 & loc_on_MAPT <= 243 ~ "Proline-rich-2",
                # loc_on_MAPT >= 244 & loc_on_MAPT <= 368 ~ "4R",
                # loc_on_MAPT >= 369 & loc_on_MAPT <= 441 ~ "C-terminal"
                loc_on_MAPT >= 1 & loc_on_MAPT <= 150 ~ "N-terminal",
                loc_on_MAPT >= 151 & loc_on_MAPT <= 243 ~ "PRR",
                loc_on_MAPT >= 244 & loc_on_MAPT <= 368 ~ "MTBR",
                loc_on_MAPT >= 369 & loc_on_MAPT <= 441 ~ "C-terminal"
            )
        )
    # colnames(df) <- colnames.ls[[i]][colnames.ls[[i]]$is_sig == TRUE, ]$label
    colnames(df) <- colnames.ls[[i]]$label
    cor_matrix <- cor(df, use = "pairwise.complete.obs", method = "pearson")
    # col_anno <- HeatmapAnnotation(Domain = factor(colnames.ls[[i]][colnames.ls[[i]]$is_sig == TRUE, ]$domain, levels = c("N-terminal", "1N", "2N", " ", "Proline-rich-1", "Proline-rich-2", "4R", "C-terminal"), ordered = TRUE),
    #                               col = list(Domain = c("N-terminal" = "grey", "1N" = "blue", "2N" = "green", " " = "grey",
    #                                                    "Proline-rich-1" = "purple", "Proline-rich-2" = "orange",
    #                                                    "4R" = "pink", "C-terminal" = "grey")),
    used_levels <- intersect(
        c("N-terminal", "PRR", "MTBR", "C-terminal"),
        # c("N-terminal", "1N", "2N", "Linker (2N-PRR)", "Proline-rich-1", "Proline-rich-2", "4R", "C-terminal"),
        unique(colnames.ls[[i]]$domain)
    )

    domain_factor <- factor(
        colnames.ls[[i]]$domain,
        levels = used_levels,
        ordered = TRUE
    )
    mycols5 <- brewer.pal(8, "Set3")

    domain_colors <- c(
        "N-terminal" = mycols5[1],
        "PRR" = mycols5[2],
        "MTBR" = mycols5[3],
        "C-terminal" = mycols5[4]
        # "N-terminal" = "grey",
        # "1N" = "#BEB8DC",
        # "2N" = "#8ECFC9",
        # "Linker (2N-PRR)" = "lightgrey",
        # "Proline-rich-1" = "#F8AC8C",
        # "Proline-rich-2" = "#9AC9DB",
        # "4R" = "#FF8884",
        # "C-terminal" = "#999999"
    )

    used_colors <- domain_colors[used_levels]

    col_anno <- HeatmapAnnotation(
        Domain = domain_factor,
        col = list(Domain = used_colors),
        show_annotation_name = FALSE,
        annotation_height = unit(5, "mm"),
        simple_anno_size_adjust = TRUE,
        gp = gpar(col = "black", lwd = 1),
        annotation_legend_param = list(
            Domain = list(
                title = "Domain",
                at = levels(domain_factor),
                labels = levels(domain_factor)
            )
        ),
        show_legend = F
    )

    row_anno <- rowAnnotation(
        Domain = domain_factor,
        col = list(Domain = used_colors),
        show_annotation_name = FALSE,
        annotation_height = unit(5, "mm"),
        simple_anno_size_adjust = TRUE,
        gp = gpar(col = "black", lwd = 1),
        annotation_legend_param = list(
            Domain = list(
                title = "Domain",
                at = levels(domain_factor),
                labels = levels(domain_factor)
            )
        )
    )

    htmap.ls[[i]] <- Heatmap(
        cor_matrix,
        name = paste0("Cor of ", names(df.ADNC.ls)[i]),
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_column_names = TRUE,
        show_row_names = TRUE,
        top_annotation = col_anno,
        left_annotation = row_anno,
        column_names_rot = 45
    )
}
pdf(
    "../Results/Tau_APP_heatmap/tau441_sig_sites_correlation_heatmap.pdf",
    width = 15,
    height = 13
)
for (i in seq_along(htmap.ls)) {
    draw(htmap.ls[[i]], padding = unit(c(10, 10, 10, 10), "mm"))
}
dev.off()
