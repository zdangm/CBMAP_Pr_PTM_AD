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
    library(clusterProfiler)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(pheatmap)
    library(impute)
})
rm(list = ls())
NUM_COV <- 56
setwd("/share/home/lik/Scripts/likun/Random_tasks/Scripts/")
RES_SUB_DIR <- "../Results/Cell_2020_heatmap/"
source(
    "/share/home/lik/Scripts/likun/Phospho/project/DEA/Scripts/helper_func.R"
)

log2_raw.729.dirs.ls <- list(
    phospho = "/data/projects/China_Brain_MultiOmics/humanBrain_Phospho/PhosPho_maintained_by_LK/DEA/Results/version_after_2025_05_23/log_impute_scale_INT/",
    ubiq = "/data/projects/China_Brain_MultiOmics/humanBrain_Ubiquitylation/Ubiquity_maintained_by_LK/DEA/Results/version_after_2025_05_23/log_impute_scale_INT/",
    ace = "/data/projects/China_Brain_MultiOmics/humanBrain_Acetylation/Acetylation_maintained_by_LK/DEA/Results/version_after_2025_06_30/log_impute_scale_INT/"
)

binary_show <- F

qs_readm(paste0(
    log2_raw.729.dirs.ls[["ace"]],
    "mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
))
ace <- log1p_raw.all.info.ls
sample_ord <- ace$intensity$id


ace <- lapply(ace, function(x) {
    x <- x[match(sample_ord, x$id), ]
    if (binary_show) {
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
        x <- is.na(x)
    } else {
        # x <- myimputeknn(x, NUM_COV)
        # x <- mymedian(x, NUM_COV)
        # x <- myINT(x, NUM_COV)
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
    }

    colnames(x) <- paste("ace", colnames(x), sep = ".")
    x
})

qs_readm(paste0(
    log2_raw.729.dirs.ls[["phospho"]],
    "mis_log2_impt_sc_rint.729.all.info.ls.QS2"
))
phospho <- log1p_raw.all.info.ls
phospho <- lapply(phospho, function(x) {
    x <- x[match(sample_ord, x$id), ]
    if (binary_show) {
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
        x <- is.na(x)
    } else {
        # x <- myimputeknn(x, NUM_COV)
        # x <- mymedian(x, NUM_COV)
        # x <- myINT(x, NUM_COV)
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
    }

    colnames(x) <- paste("phospho", colnames(x), sep = ".")
    x
})

qs_readm(paste0(
    log2_raw.729.dirs.ls[["ubiq"]],
    "mis_log2_impt_sc_rint.729.all.info.ls.QS2"
))
ubiq <- log1p_raw.all.info.ls
ubiq <- lapply(ubiq, function(x) {
    x <- x[match(sample_ord, x$id), ]

    if (binary_show) {
        x[, seq.int(1, ncol(x) - NUM_COV)] <- is.na(x[, seq.int(
            1,
            ncol(x) - NUM_COV
        )])
    } else {
        # x <- myimputeknn(x, NUM_COV)
        # x <- mymedian(x, NUM_COV)
        # x <- myINT(x, NUM_COV)
    }

    colnames(x)[seq.int(1, ncol(x) - NUM_COV)] <- paste(
        "ubiq",
        colnames(x)[seq.int(1, ncol(x) - NUM_COV)],
        sep = "."
    )
    x
})


log1p_raw.all.info.ls <- mapply(
    function(a, p, u) {
        x <- cbind(a, p, u)
        x
    },
    ace,
    phospho,
    ubiq,
    SIMPLIFY = F
)


df.hc.ls <- lapply(log1p_raw.all.info.ls, function(x) {
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


#------For ADNC catagory-------#
#------For ADNC catagory-------#
#------For ADNC catagory-------#
#------For ADNC catagory-------#
df.ADNC.case.ls <- lapply(log1p_raw.all.info.ls, function(x) {
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
    x <- x %>% dplyr::filter(!(ADNC_LMH == 0 & A_beta_0_3 == 1))
    x
})
lapply(df.ADNC.ls, dim)
df.ADNC.info <- df.ADNC.ls[[1]][, -seq.int(1, ncol(df.ADNC.ls[[1]]) - NUM_COV)]


df.ADNC <- df.ADNC.ls[[1]]


phospho.sig <- read.xlsx(paste0(
    log2_raw.729.dirs.ls[["phospho"]],
    "All.sites.All.traits.intensity.limma.xlsx"
))
phospho.sig <- phospho.sig[phospho.sig$ADNC.LMH.num.adj.P.Val < 0.05, ]$id
phospho.sig <- paste("phospho", phospho.sig, sep = ".")
length(phospho.sig)
ubiq.sig <- read.xlsx(paste0(
    log2_raw.729.dirs.ls[["ubiq"]],
    "All.sites.All.traits.intensity.limma.xlsx"
))
ubiq.sig <- ubiq.sig[ubiq.sig$ADNC.LMH.num.adj.P.Val < 0.05, ]$id
ubiq.sig <- paste("ubiq", ubiq.sig, sep = ".")
length(ubiq.sig)
ace.sig <- read.xlsx(paste0(
    log2_raw.729.dirs.ls[["ace"]],
    "All.sites.All.traits.intensity.limma.xlsx"
))
ace.sig <- ace.sig[ace.sig$ADNC.LMH.num.adj.P.Val < 0.05, ]$id
ace.sig <- paste("ace", ace.sig, sep = ".")
length(ace.sig)
head(ace.sig)

# df.ADNC <- df.ADNC %>%
#     dplyr::select(any_of(c(phospho.sig, ubiq.sig, ace.sig)))
dim(df.ADNC)
dim(df.ADNC.info)

#  OR
df.ADNC <- df.ADNC %>%
    dplyr::select(contains("MAPT"))
dim(df.ADNC)
dim(df.ADNC.info)
domain.order <- str_split_i(colnames(df.ADNC), "\\.", 3)
df.ADNC <- df.ADNC[, gtools::mixedorder(domain.order)]
df.ADNC.info <- df.ADNC.info %>% arrange(ADNC_LMH)
df.ADNC <- df.ADNC[match(df.ADNC.info$id, rownames(df.ADNC)), ]
sum(rownames(df.ADNC) == df.ADNC.info$id)

palettes <- c("Blues", "Greens", "Oranges", "Purples")

#
colors_list <- lapply(palettes, function(p) brewer.pal(4, p))

mycols <- brewer.pal(3, "Pastel1")
mycols5 <- brewer.pal(8, "Set3")


df.ADNC <- df.ADNC[,
    sapply(
        str_split_i(colnames(df.ADNC), "_", 2),
        from_longest_isof_to_clinical_isof
    ) !=
        "not found on MAPT-8"
]
sites <- str_split_i(
    sapply(
        str_split_i(colnames(df.ADNC), "_", 2),
        from_longest_isof_to_clinical_isof
    ),
    "_",
    2
)
sites <- as.numeric(sites)
dim(df.ADNC)
domain <- cut(
    sites,
    breaks = c(0, 150, 243, 368, 441),
    labels = c("N-terminal", "PRR", "MTBR", "C-terminal"),
    right = TRUE
)
# domain <- cut(
#             sites,
#             breaks = c(0, 44, 73, 102, 150, 197, 243, 368, 441),
#             labels = c("N-terminal", "1N", "2N", "spacer", "Proline-rich-1", "Proline-rich-2", "4R", "C-terminal"),
#             right = TRUE
#             )

ha_row <- HeatmapAnnotation(
    PTM = str_split_i(colnames(df.ADNC), "\\.", 1),
    Domain = domain,
    col = list(
        PTM = c("phospho" = mycols[1], "ubiq" = mycols[2], "ace" = mycols[3]),
        Domain = c(
            "N-terminal" = mycols5[1],
            "PRR" = mycols5[2],
            "MTBR" = mycols5[3],
            "C-terminal" = mycols5[4]
        )
        # domain = c("N-terminal" = mycols5[1], "1N" = mycols5[2], "2N" = mycols5[3], "spacer" = mycols5[4], "Proline-rich-1" = mycols5[5], "Proline-rich-2" = mycols5[6], "4R" = mycols5[7], "C-terminal" = mycols5[8])
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 16),
    annotation_height = unit(rep(8, 2), "mm"),
    simple_anno_size_adjust = TRUE,
    gp = gpar(col = "darkgrey", lwd = .5),
    annotation_legend_param = list(
        title_gp = gpar(fontsize = 16),
        labels_gp = gpar(fontsize = 14)
    )
)

mycolors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(4)

ha_col <- rowAnnotation(
    # Sex_male = df.ADNC.info$sex_male,
    # bank = df.ADNC.info$bank,
    # Dementia = df.ADNC.info$diag_dementia,
    A_beta = df.ADNC.info$A_beta_0_3,
    B_nft = df.ADNC.info$B_NFT_braak_0_3,
    # C_cerad = df.ADNC.info$C_cerad_0_3,
    ADNC_LMH = df.ADNC.info$ADNC_LMH,
    # ADNC = df.ADNC.info$ADNC,
    # digest_batch = df.ADNC.info$digest_batch,
    # CVD = df.ADNC.info$CVD,
    # LBD = df.ADNC.info$LBD,
    # LATE = df.ADNC.info$LATE,
    # ARTAG = df.ADNC.info$ARTAG,
    col = list(
        ADNC_LMH = c(
            "0" = mycolors[1],
            "1" = mycolors[2],
            "2" = mycolors[3],
            "3" = mycolors[4]
        ),
        B_nft = c(
            "0" = colors_list[[2]][1],
            "1" = colors_list[[2]][2],
            "2" = colors_list[[2]][3],
            "3" = colors_list[[2]][4]
        ),
        # C_cerad = c("0" = colors_list[[3]][1], "1" = colors_list[[3]][2], "2" = colors_list[[3]][3], "3" = colors_list[[3]][4]),
        A_beta = c(
            "0" = colors_list[[4]][1],
            "1" = colors_list[[4]][2],
            "2" = colors_list[[4]][3],
            "3" = colors_list[[4]][4]
        )
        # ADNC = c("0" = "white", "1" = "black"),
        # Sex_male = c("0" = "white", "1" = "black"),
        # bank = c("csu" = "blue", "pumc" = "white", "zju" = "red"),
        # Dementia = c("0" = "white", "1" = "black")
        # CVD = c("0" = "white", "1" = "black"),
        # LBD = c("0" = "white", "1" = "black"),
        # LATE = c("0" = "white", "1" = "black"),
        # ARTAG = c("0" = "white", "1" = "black"),
        # cluster = c("1" = "blue", "2" = "red", "3" = "green"),
        # digest_batch = c("A1" = "blue", "A2" = "red", "A3" = "green", "A4" = "black")
    ),
    show_annotation_name = F,
    annotation_name_gp = gpar(fontsize = 16),
    annotation_width = unit(rep(8, 3), "mm"),
    simple_anno_size_adjust = T,
    # width = unit(39, "mm"),

    # simple_anno_size = unit(15, "mm"),
    annotation_name_rot = 45,
    annotation_legend_param = list(
        title_gp = gpar(fontsize = 16),
        labels_gp = gpar(fontsize = 14)
    )
)
getwd()


if (binary_show) {
    rn <- rownames(df.ADNC)
    df.ADNC <- as.data.frame(lapply(df.ADNC, function(x) {
        x <- if_else(x, 0, 1)
        # x <- factor(x, levels = c(0, 1), labels = c("identified", "not identified"))
        x
    }))
    rownames(df.ADNC) <- rn
}
pdf(
    sprintf("../Results/Cell_2020_heatmap/MAPT_all_sites%s.pdf", "no_cluster"),
    width = 15,
    height = 15
)
ht <- Heatmap(
    as.matrix(df.ADNC),
    name = "TAU",
    col = circlize::colorRamp2(
        c(min(df.ADNC), max(df.ADNC)),
        c("white", "red")
    ),

    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    # clustering_distance_rows = "manhattan",
    # clustering_distance_columns = "manhattan",
    # clustering_method_rows = "ward.D2",
    # clustering_method_columns = "ward.D2",
    show_column_names = FALSE,
    show_row_names = FALSE,
    left_annotation = ha_col,
    top_annotation = ha_row,
    heatmap_height = unit(8, "inches"),
    heatmap_width = unit(8, "inches"),
    # row_dend_width = unit(1, "inches"),
    heatmap_legend_param = list(
        title_gp = gpar(fontsize = 0),
        labels_gp = gpar(fontsize = 14)
    ),
    # width = unit(8, "inches"),
    # height = unit(8, "inches"),
    border = FALSE
)
draw(ht)
dev.off()
