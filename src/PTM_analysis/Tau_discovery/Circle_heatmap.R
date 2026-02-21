#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

suppressPackageStartupMessages({
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
    library(tidyverse)
})
rm(list = ls())
setwd("/share/home/lik/Scripts/likun/Random_tasks/Scripts/")
source(
    "/share/home/lik/Scripts/likun/Phospho/project/DEA/Scripts/helper_func.R"
)
conflicted::conflict_prefer_all("dplyr")
RES_SUB_DIR.ls <- list(
    phospho = "/data/projects/China_Brain_MultiOmics/humanBrain_Phospho/PhosPho_maintained_by_LK/DEA/Results/version_after_2025_05_23/",
    ubiq = "/data/projects/China_Brain_MultiOmics/humanBrain_Ubiquitylation/Ubiquity_maintained_by_LK/DEA/Results/version_after_2025_05_23/",
    ace = "/data/projects/China_Brain_MultiOmics/humanBrain_Acetylation/Acetylation_maintained_by_LK/DEA/Results/version_after_2025_06_30/"
)
NUM_COV <- 56

is_tau <- T

#--------------------------------------------------Begin---------------------------------------------------------#
#   这一部分是用missing rate < 0.5的位点来画图,impute 且 纵向scale int的intensity来画热图，得到MAPT.ls.1
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
    for (j in seq_along(log1p_raw.all.info.ls)) {
        x <- log1p_raw.all.info.ls[[j]]
        # x <- myimputeknn(x, NUM_COV =  56); x <- mymedian(x, NUM_COV = 56); x <- myINT(x, NUM_COV = 56)
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
lapply(MAPT.ls.1, dim)

#------------------------------------------------------END-------------------------------------------------------#

#---------------------------------#
#-------prepare covariates--------#
#---------------------------------#
batch_info <- read.xlsx(
    "/share/home/lik/Scripts/likun/sample_meta_info/【样本明细】966标+21标MIX（共987标）人脑组织4D fastDIA蛋白组、磷酸化、乙酰化、泛素化项目正式样本明细20241121.xlsx",
    sheet = 3
)

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
dim(covariates)
# 729 54
rm(sample_info, batch_info)
gc()


#----------------------------------------------------------#
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
df.ADNC.ls <- mapply(
    function(x, y) {
        rbind(x, y)
    },
    df.hc.ls,
    df.ADNC.case.ls,
    SIMPLIFY = F
)


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
df.ADNC.avg.ls <- list()
i <- 1
for (i in seq_along(df.ADNC.ls)) {
    df <- df.ADNC.ls[[i]]
    df <- df %>%
        group_by(ADNC_LMH) %>%
        summarise(across(matches("^phospho|^ace|^ubiq"), median, na.rm = TRUE))
    df.ADNC.avg.ls[[names(df.ADNC.ls)[i]]] <- df
}
lapply(df.ADNC.avg.ls, dim)
head(df.ADNC.avg.ls$intensity)
i <- 1
j <- 1

for (i in seq_along(df.ADNC.avg.ls)) {
    df <- df.ADNC.avg.ls[[i]]
    df <- as.data.frame(df, check.names = F)
    df <- df %>% column_to_rownames("ADNC_LMH")
    df <- t(df)
    df <- as.data.frame(df, check.names = F)
    df <- df %>%
        rownames_to_column("identifier") %>%
        separate(
            identifier,
            into = c("ptm", "gene_name", "identifier"),
            sep = "\\."
        ) %>%
        unite(identifier, gene_name, identifier, sep = ".")
    df.ls <- df %>% group_by(ptm) %>% group_split(.keep = TRUE)
    df.ptm.ls <- list()
    for (j in seq_along(df.ls)) {
        df.ptm <- df.ls[[j]]

        df.ptm <- df.ptm %>% column_to_rownames("identifier")
        df.ptm <- get_gene_name(df.ptm, PTM = df.ptm$ptm[1])
        if (is_tau) {
            df.ptm$loc_on_tau441 <- sapply(
                df.ptm$for_volcano,
                from_longest_isof_to_clinical_isof
            )
            df.ptm <- df.ptm %>%
                filter(
                    loc_on_tau441 != "not found on MAPT-8",
                    !is.na(loc_on_tau441)
                )
            df.ptm$label <- paste0(
                str_split_i(df.ptm$loc_on_tau441, "_", 3),
                str_split_i(df.ptm$loc_on_tau441, "_", 2)
            )
            df.ptm$loc_on_tau441 <- str_split_i(df.ptm$loc_on_tau441, "_", 2)
        } else {
            df.ptm <- df.ptm %>%
                separate(
                    for_volcano,
                    into = c("Gene_name", "loc_on_tau441", "aa"),
                    sep = "_"
                )
            df.ptm$label <- paste0(df.ptm$aa, df.ptm$loc_on_tau441)
        }
        df.ptm <- df.ptm %>%
            select(ptm, HC, Low, Moderate, High, label, loc_on_tau441)
        rownames(df.ptm) <- NULL
        df.ptm.ls[[j]] <- df.ptm
    }
    df <- bind_rows(df.ptm.ls)
    df.ADNC.avg.ls[[i]] <- df
}
for (i in seq_along(df.ADNC.avg.ls)) {
    df <- df.ADNC.avg.ls[[i]]
    boolvec <- rowSums(is.na(df)) < 2
    message(any(boolvec))
    df <- df[boolvec, ]
    df.ADNC.avg.ls[[i]] <- df
}
lapply(df.ADNC.avg.ls, dim)
for (i in seq_along(df.ADNC.avg.ls)) {
    df <- df.ADNC.avg.ls[[i]]
    df$loc_on_tau441 <- as.numeric(df$loc_on_tau441)
    if (is_tau) {
        # df$domain <- cut(
        #     df$loc_on_tau441,
        #     breaks = c(0, 44, 73, 102, 150, 197, 243, 368, 441),
        #     labels = c("N-terminal", "1N", "2N", " ", "Proline-rich-1", "Proline-rich-2", "4R", "C-terminal"),
        #     right = TRUE
        #     )
        df$domain <- cut(
            df$loc_on_tau441,
            breaks = c(0, 150, 243, 368, 441),
            labels = c("N-terminal", "PRR", "MTBR", "C-terminal"),
            right = TRUE
        )
    } else {
        df$domain <- cut(
            df$loc_on_tau441,
            breaks = c(0, 671, 713, 770),
            labels = c("APPs-beta", "A-beta", "AICD"),
            right = TRUE
        )
    }
    df.ADNC.avg.ls[[i]] <- df
}
for (i in seq_along(df.ADNC.avg.ls)) {
    df <- df.ADNC.avg.ls[[i]]
    val_cols <- df %>% select(HC:High) %>% as.matrix()
    val_cols <- t(apply(val_cols, 1, scale))
    val_cols <- as.data.frame(val_cols, check.names = F)
    colnames(val_cols) <- c("HC", "Low", "Moderate", "High")
    df <- df %>% select(-HC, -Low, -Moderate, -High)
    df <- cbind(df, val_cols)
    df.ADNC.avg.ls[[i]] <- df
}
domain.ls <- list()
for (i in seq_along(df.ADNC.avg.ls)) {
    df <- df.ADNC.avg.ls[[i]]
    df <- df %>% arrange(loc_on_tau441)
    if (is_tau) {
        domain.ls[[i]] <- factor(
            df$domain,
            levels = c("N-terminal", "PRR", "MTBR", "C-terminal"),
            ordered = TRUE
        )
        # domain.ls[[i]] <- factor(df$domain, levels = c("N-terminal", "1N", "2N", " ", "Proline-rich-1", "Proline-rich-2", "4R", "C-terminal"), ordered = TRUE)
    } else {
        domain.ls[[i]] <- factor(
            df$domain,
            levels = c("APPs-beta", "A-beta", "AICD"),
            ordered = TRUE
        )
    }
    domain.ls[[i]] <- droplevels(domain.ls[[i]])
    rownames(df) <- paste(df$ptm, df$label, sep = ".")
    df <- df %>% select(-ptm, -loc_on_tau441, -label, -domain)
    df <- as.matrix(df)
    df.ADNC.avg.ls[[i]] <- df
}

col_color <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) #
missing_lt0.5.ls <- list(intensity = NULL, no_intensity = NULL)
sig.sites.ls <- list(intensity = NULL, no_intensity = NULL)
All.sites.df.ls <- list(intensity = NULL, no_intensity = NULL)
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
        # all_sites.dir <- file.path(all_sites.dir, "version_after_2025_05_23/log_impute_scale_INT/All.sites.limma.QS2")
        qs_readm(all_sites.dir)
        site_names <- lapply(All.sites, function(x) {
            x <- rownames(x)
            if (is_tau) {
                x <- x[str_detect(x, "^MAPT\\.P10636")]
                x <- data.frame(id = x, check.names = F)
                x <- x %>%
                    separate(
                        id,
                        into = c("pr", "pos", "aa"),
                        sep = "_",
                        remove = F
                    )
                x$pos <- sapply(x$pos, from_longest_isof_to_clinical_isof)
                x <- x %>%
                    filter(pos != "not found on MAPT-8") %>%
                    separate(
                        pos,
                        into = c("pos1", "pos"),
                        sep = "_",
                        remove = F
                    )
                x$label <- paste0(names(RES_SUB_DIR.ls)[j], ".", x$aa, x$pos)
            } else {
                x <- x[str_detect(x, "^APP\\.P05067")]
                x <- data.frame(id = x, check.names = F)
                x <- x %>%
                    separate(
                        id,
                        into = c("pr", "pos", "aa"),
                        sep = "_",
                        remove = F
                    )
                x$label <- paste0(names(RES_SUB_DIR.ls)[j], ".", x$aa, x$pos)
            }
            x$label
        })
        allsig.sites <- lapply(All.sites, function(x) {
            if (!is_tau) {
                x <- x %>%
                    filter(ADNC.LMH.num.P.Value < 0.05) %>%
                    select(
                        ADNC.LMH.num.P.Value,
                        braak.num.P.Value,
                        ADNC.LMH.num.adj.P.Val
                    ) %>%
                    rownames_to_column("identifier") %>%
                    filter(str_detect(identifier, "^APP\\.P05067"))
            } else {
                x <- x %>%
                    filter(ADNC.LMH.num.adj.P.Val < 0.05) %>% #
                    select(
                        ADNC.LMH.num.P.Value,
                        braak.num.P.Value,
                        ADNC.LMH.num.adj.P.Val
                    ) %>%
                    rownames_to_column("identifier") %>%
                    filter(str_detect(identifier, "^MAPT\\.P10636"))
            }
            if (nrow(x) == 0) {
                message(
                    "No significant sites found in ",
                    names(RES_SUB_DIR.ls)[j],
                    " for ",
                    names(sig.sites.ls)[i]
                )
                return(NULL)
            }
            x$identifier <- paste(
                names(RES_SUB_DIR.ls)[j],
                x$identifier,
                sep = "."
            )
            x
        })
        All.sites <- lapply(All.sites, function(x) {
            if (!is_tau) {
                x <- x %>%
                    select(
                        ADNC.LMH.num.P.Value,
                        braak.num.P.Value,
                        ADNC.LMH.num.adj.P.Val
                    ) %>%
                    rownames_to_column("identifier") %>%
                    filter(str_detect(identifier, "^APP\\.P05067"))
            } else {
                x <- x %>%
                    select(
                        ADNC.LMH.num.P.Value,
                        braak.num.P.Value,
                        ADNC.LMH.num.adj.P.Val
                    ) %>%
                    rownames_to_column("identifier") %>%
                    filter(str_detect(identifier, "^MAPT\\.P10636"))
            }
            x$identifier <- paste(
                names(RES_SUB_DIR.ls)[j],
                x$identifier,
                sep = "."
            )
            x
        })
        All.sites.df.ls[[i]] <- rbind(
            All.sites.df.ls[[i]],
            All.sites[[i]][,
                c("identifier", "ADNC.LMH.num.adj.P.Val"),
                drop = FALSE
            ]
        )
        site_names <- site_names[[i]]
        missing_lt0.5.ls[[i]] <- c(missing_lt0.5.ls[[i]], site_names)
        if (is.null(All.sites)) {
            next
        }
        sig.sites.ls[[i]] <- rbind(
            sig.sites.ls[[i]],
            allsig.sites[[i]][,
                c("identifier", "ADNC.LMH.num.adj.P.Val"),
                drop = FALSE
            ]
        )
    }
}
for (i in seq_along(sig.sites.ls)) {
    df <- sig.sites.ls[[i]]

    df$ptm <- str_split_i(df$identifier, "\\.", 1)
    df$pos <- str_split_i(df$identifier, "_", 2)
    if (is_tau) {
        df$pos <- sapply(df$pos, from_longest_isof_to_clinical_isof)
        df <- df %>% filter(pos != "not found on MAPT-8")
        df$pos <- str_split_i(df$pos, "_", 2)
    }
    df$aa <- str_split_i(df$identifier, "_", 3)
    df$label <- paste0(df$ptm, ".", df$aa, df$pos)
    sig.sites.ls[[i]] <- df

    df <- All.sites.df.ls[[i]]

    df$ptm <- str_split_i(df$identifier, "\\.", 1)
    df$pos <- str_split_i(df$identifier, "_", 2)
    if (is_tau) {
        df$pos <- sapply(df$pos, from_longest_isof_to_clinical_isof)
        df <- df %>% filter(pos != "not found on MAPT-8")
        df$pos <- str_split_i(df$pos, "_", 2)
    }
    df$aa <- str_split_i(df$identifier, "_", 3)
    df$label <- paste0(df$ptm, ".", df$aa, df$pos)
    All.sites.df.ls[[i]] <- df
}

sig.col.ls <- list(intensity = NULL, no_intensity = NULL)
for (i in seq_along(sig.col.ls)) {
    rn <- rownames(df.ADNC.avg.ls[[i]])
    bv <- rn %in% sig.sites.ls[[i]]$label
    sig.col.ls[[i]] <- if_else(
        bv,
        "red",
        if_else(rn %in% missing_lt0.5.ls[[i]], "black", "lightgrey")
    )
}

if (is_tau) {
    qs_readm("pure_assoc_p.QS2")
    qs_readm("sequential_assoc_p.QS2")
    # pure_assoc_p <- lapply(pure_assoc_p, function(x) {
    #     x <- unlist(x)
    #     x <- x[x < 0.05]
    #     res <- data.frame(
    #         id = names(x),
    #         aa = str_split_i(names(x), "_", 3),
    #         pos = str_split_i(names(x), "_", 2),
    #         ptm = str_split_i(names(x), "\\.", 1)
    #     )
    #     res$label <- paste0(res$ptm, ".", res$aa, res$pos)
    #     res$label

    # })
    # or
    impute_p <- data.frame(label = rownames(df.ADNC.avg.ls[[1]]), p = 1)
    pure_assoc_p <- lapply(pure_assoc_p, function(x) {
        res <- data.frame(
            id = names(x),
            aa = str_split_i(names(x), "_", 3),
            pos = str_split_i(names(x), "_", 2),
            ptm = str_split_i(names(x), "\\.", 1),
            raw_p = unlist(x)
        )
        res$label <- paste0(res$ptm, ".", res$aa, res$pos)
        res <- res[, c("label", "raw_p")]
        res <- res %>% right_join(impute_p, by = "label")
        res$raw_p <- if_else(is.na(res$raw_p), 1, res$raw_p)
        # res$adj_p <- p.adjust(res$raw_p, method = "BH")
        res$adj_p <- res$raw_p
        message("using raw p ")
        res
    })
    # or
    # pure_assoc_p <- lapply(pure_assoc_p, function(x) {
    #     x <- unlist(x)
    #     res <- data.frame(
    #         id = names(x),
    #         aa = str_split_i(names(x), "_", 3),
    #         pos = str_split_i(names(x), "_", 2),
    #         ptm = str_split_i(names(x), "\\.", 1),
    #         adj_p = x
    #     )
    #     res$label <- paste0(res$ptm, ".", res$aa, res$pos)
    #     res
    # })

    adj.sig.df <- data.frame(
        id = rownames(df.ADNC.avg.ls[[1]]),
        check.names = F
    )

    # adj.sig.df[["no adj"]] <- if_else(adj.sig.df$id %in% sig.sites.ls[[1]]$label, 1, 0)
    # adj.sig.df[["adj top1"]] <- if_else(adj.sig.df$id %in% pure_assoc_p$top1, 1, 0)
    # adj.sig.df[["adj top2"]] <- if_else(adj.sig.df$id %in% pure_assoc_p$top2, 1, 0)
    # adj.sig.df[["adj top3"]] <- if_else(adj.sig.df$id %in% pure_assoc_p$top3, 1, 0)
    # adj.sig.df[["adj top4"]] <- if_else(adj.sig.df$id %in% pure_assoc_p$top4, 1, 0)
    # or
    adj.sig.df[["no adj"]] <- 1
    adj.sig.df[["adj ub.K311"]] <- 1
    adj.sig.df[["+ ub.K321"]] <- 1
    adj.sig.df[["+ ub.K254"]] <- 1
    adj.sig.df[["+ ac.K274"]] <- 1
    # adj.sig.df[["adj S262"]] <- 1
    # adj.sig.df[["+ S289"]] <- 1
    # adj.sig.df[["+ T212"]] <- 1
    # adj.sig.df[["+ T217"]] <- 1
    for (i in seq_along(adj.sig.df$id)) {
        no_adj_idx <- match(adj.sig.df$id[i], All.sites.df.ls[[1]]$label)

        adj.sig.df[i, "no adj"] <- if_else(
            is.na(no_adj_idx),
            1,
            All.sites.df.ls[[1]][no_adj_idx, "ADNC.LMH.num.adj.P.Val"]
        )
        top1_idx <- match(adj.sig.df$id[i], pure_assoc_p$top1$label)
        adj.sig.df[i, "adj ub.K311"] <- if_else(
            is.na(top1_idx),
            1,
            pure_assoc_p$top1[top1_idx, "adj_p"]
        )

        top2_idx <- match(adj.sig.df$id[i], pure_assoc_p$top2$label)
        adj.sig.df[i, "+ ub.K321"] <- if_else(
            is.na(top2_idx),
            1,
            pure_assoc_p$top2[top2_idx, "adj_p"]
        )

        top3_idx <- match(adj.sig.df$id[i], pure_assoc_p$top3$label)
        adj.sig.df[i, "+ ub.K254"] <- if_else(
            is.na(top3_idx),
            1,
            pure_assoc_p$top3[top3_idx, "adj_p"]
        )

        top4_idx <- match(adj.sig.df$id[i], pure_assoc_p$top4$label)
        adj.sig.df[i, "+ ac.K274"] <- if_else(
            is.na(top4_idx),
            1,
            pure_assoc_p$top4[top4_idx, "adj_p"]
        )
        # top1_idx <- match(adj.sig.df$id[i], pure_assoc_p$top1$label)
        # adj.sig.df[i, "adj S262"] <- if_else(is.na(top1_idx), 1,
        #     pure_assoc_p$top1[top1_idx, "adj_p"])

        # top2_idx <- match(adj.sig.df$id[i], pure_assoc_p$top2$label)
        # adj.sig.df[i, "+ S289"] <- if_else(is.na(top2_idx), 1,
        #     pure_assoc_p$top2[top2_idx, "adj_p"])

        # top3_idx <- match(adj.sig.df$id[i], pure_assoc_p$top3$label)
        # adj.sig.df[i, "+ T212"] <- if_else(is.na(top3_idx), 1,
        #     pure_assoc_p$top3[top3_idx, "adj_p"])

        # top4_idx <- match(adj.sig.df$id[i], pure_assoc_p$top4$label)
        # adj.sig.df[i, "+ T217"] <- if_else(is.na(top4_idx), 1,
        #     pure_assoc_p$top4[top4_idx, "adj_p"])
    }

    adj.sig.df <- adj.sig.df %>%
        column_to_rownames("id")
    head(adj.sig.df)
}


source("./circos.heatmap2.R")
sig_mat <- matrix(
    ifelse(as.matrix(adj.sig.df) < 0.05, "*", ""),
    nrow = nrow(adj.sig.df),
    ncol = ncol(adj.sig.df)
)

rownames(sig_mat) <- rownames(adj.sig.df)
sig_mat <- sig_mat[, seq.int(ncol(sig_mat), 1)] # reverse the order of columns
colnames(sig_mat) <- colnames(adj.sig.df)
adj.sig.df <- -log10(as.matrix(adj.sig.df))
maxp <- max(adj.sig.df, na.rm = TRUE)

conflicted::conflicts_prefer(base::unname)


pdf(
    "../Results/Tau_APP_heatmap/tau441_circle_heatmap_imputed_scaled_INT_adj_top_sites.pdf",
    width = 10,
    height = 10
)
circos.par(
    start.degree = 65,
    gap.after = c(rep(1, length(levels(domain.ls[[1]])) - 1), 25)
) #让分裂的一个口大一点，可以添加行信息

circos.heatmap(
    adj.sig.df,
    split = domain.ls[[1]],
    col = colorRamp2(
        c(maxp, 0),
        str_split_fixed("red, white", ", ", n = Inf)[1, ]
    ),
    na.col = "darkgrey",
    rownames.cex = 0.9,
    dend.side = "none",
    rownames.side = "outside",
    cluster = FALSE,
    cell.border = NA,
    rownames.col = "black",
    show.sector.labels = F,
    bg.border = "black",
    mat.sig = sig_mat
)


circlize::circos.track(
    track.index = 2,
    panel.fun = function(x, y) {
        if (CELL_META$sector.numeric.index == length(levels(domain.ls[[1]]))) {
            # the last sector
            cn = colnames(adj.sig.df)
            n = length(cn)
            y_pos = c(5, 4, 3, 2, 1) - 0.5
            circos.text(
                rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                y_pos,
                cn,
                cex = 1,
                adj = c(0, 0.5),
                facing = "inside"
            )
        }
    },
    bg.border = NA
)
circlize::circos.heatmap(
    df.ADNC.avg.ls[[1]],
    split = domain.ls[[1]],
    col = col_color,
    na.col = "darkgrey",
    # rownames.cex = ,
    dend.side = "none",
    rownames.side = "none",
    cluster = FALSE,
    # rownames.col = sig.col.ls[[1]],
    show.sector.labels = F,
    bg.border = "black"
)

circlize::circos.track(
    track.index = 3,
    panel.fun = function(x, y) {
        if (CELL_META$sector.numeric.index == length(levels(domain.ls[[1]]))) {
            # the last sector
            cn = colnames(df.ADNC.avg.ls[[1]])
            n = length(cn)
            y_pos = c(4, 3, 2, 1) - 0.5
            circos.text(
                rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                y_pos,
                cn,
                cex = 1,
                adj = c(0, 0.5),
                facing = "inside"
            )
        }
    },
    bg.border = NA
)


circlize::circos.track(
    track.index = 2,
    panel.fun = function(x, y) {
        circos.text(
            CELL_META$xcenter,
            CELL_META$cell.ylim[-25] +
                convert_y(-30, "mm"),
            CELL_META$sector.index,
            facing = "bending.outside",
            cex = .9,
            adj = c(0.5, 0),
            niceFacing = T
        )
    },
    bg.border = NA
)


lg1 <- Legend(
    title = "Intensity",
    col_fun = col_color,
    direction = "vertical",
    title_position = "topleft"
)
lg2 <- Legend(
    title = bquote("-log"[10] * "P"),
    col_fun = colorRamp2(
        c(maxp, 0),
        str_split_fixed("red, white", ", ", n = Inf)[1, ]
    ),
    direction = "vertical",
    title_position = "topleft"
)

packed_legend <- packLegend(lg1, lg2, direction = "vertical")

draw(
    packed_legend,
    x = unit(1, "npc") - unit(2, "mm"),
    y = unit(4, "mm"),
    just = c("right", "bottom")
)

circos.clear()
dev.off()
