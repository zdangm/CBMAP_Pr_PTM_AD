suppressPackageStartupMessages({
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
    library(impute)
    library(ggrepel)
    library(clusterProfiler)
    library(tidyverse)
    library(conflicted)
    options(tibble.width = Inf)
})
rm(list = ls())
conflict_prefer_all("dplyr")
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
RES_SUB_DIR <- paste0("../Results/", PTMs[i])
NUM_COV <- 56


#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------helper functions-------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
source(
    "helper_func.R"
)
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------DONE--------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

ppline <- "log_impute_scale_INT"


RES_SUB_DIR <- paste0(RES_SUB_DIR, ppline)
RES_SUB_DIR <- paste0(RES_SUB_DIR, "/")
dir.exists(RES_SUB_DIR)

if (dir.exists(RES_SUB_DIR)) {
    str <- sprintf("Directory %s created successfully", RES_SUB_DIR)
    message(str)
}
if (PTMs[i] == "ace") {
    qs_readm(paste0(
        RES_SUB_DIR,
        "mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
    ))
} else {
    qs_readm(paste0(RES_SUB_DIR, "mis_log2_impt_sc_rint.729.all.info.ls.QS2"))
}
lapply(log1p_raw.all.info.ls, dim)

#----------------------------------------------------For ADNC---------------------------------------------------#
#----------------------------------------------------For ADNC---------------------------------------------------#
#----------------------------------------------------For ADNC---------------------------------------------------#
df.hc.ls <- lapply(log1p_raw.all.info.ls, function(x) {
    df.hc <- x %>%
        dplyr::filter(
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

df.ADNC.case.ls <- lapply(log1p_raw.all.info.ls, function(x) {
    df.ADNC.case <- x %>%
        dplyr::filter(
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
df.ADNC.ls <- mapply(
    function(x, y) {
        res <- rbind(x, y)
    },
    df.hc.ls,
    df.ADNC.case.ls,
    SIMPLIFY = F
)
lapply(df.ADNC.ls, dim)
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------#

#----------------------------------------------For braak NFT stage----------------------------------------------#
#----------------------------------------------For braak NFT stage----------------------------------------------#
#----------------------------------------------For braak NFT stage----------------------------------------------#
df.braak.ls <- lapply(log1p_raw.all.info.ls, function(x) {
    df.braak <- x %>%
        dplyr::filter(
            OTHER == 0,
            diag_ALS == 0,
            diag_SCZ == 0,
            diag_epilepsy == 0,
            !is.na(Braak_NFT_stage),
            diag_others == 0
        )
    df.braak
})
lapply(df.braak.ls, dim)

table(df.braak.ls[[1]]$Braak_NFT_stage, useNA = "always")
table(df.braak.ls[[1]]$B_NFT_braak_0_3, useNA = "always")


#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------AB scores---------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
df.a.score.ls <- lapply(log1p_raw.all.info.ls, function(x) {
    df.braak <- x %>%
        dplyr::filter(
            OTHER == 0,
            diag_ALS == 0,
            diag_SCZ == 0,
            diag_epilepsy == 0,
            !is.na(A_beta_0_3), ## 2025-05-23  maybe changed
            diag_others == 0

            # bank != "csu",
        )
    df.braak
})
lapply(df.a.score.ls, dim)

df.b.score.ls <- lapply(log1p_raw.all.info.ls, function(x) {
    df.braak <- x %>%
        dplyr::filter(
            OTHER == 0,
            diag_ALS == 0,
            diag_SCZ == 0,
            diag_epilepsy == 0,
            !is.na(B_NFT_braak_0_3),
            diag_others == 0
        )
    df.braak
})
lapply(df.b.score.ls, dim)


#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------DONE-------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------Diff Analysis--------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
cov_names <- c("sex_male", "age", "PMD", "RIN") # include IP_batch only when sensitivity analysis is needed

DEA.a.score.limma.res.numeric.ls <- diff_analysis(
    c("A_beta_0_3", cov_names),
    df.c.score.ls,
    method = "limma",
    var_type = "numeric",
    NUM_COV = NUM_COV
)
DEA.b.score.limma.res.numeric.ls <- diff_analysis(
    c("B_NFT_braak_0_3", cov_names),
    df.b.score.ls,
    method = "limma",
    var_type = "numeric",
    NUM_COV = NUM_COV
)
DEA.braak.limma.res.numeric.ls <- diff_analysis(
    c("Braak_NFT_stage", cov_names),
    df.braak.ls,
    method = "limma",
    var_type = "numeric",
    NUM_COV = NUM_COV
)
DEA.ADNC.LMH.limma.res.numeric.ls <- diff_analysis(
    c("ADNC_LMH", cov_names),
    df.ADNC.ls,
    method = "limma",
    var_type = "numeric",
    NUM_COV = NUM_COV
)
head(tmp)
qs_savem(
    DEA.a.score.limma.res.numeric.ls,
    DEA.b.score.limma.res.numeric.ls,
    DEA.ADNC.LMH.limma.res.numeric.ls,
    DEA.braak.limma.res.numeric.ls,
    file = paste0(RES_SUB_DIR, "DEA.limma.res.ls.QS2")
)


df.adj.ADNC.ls <- lapply(df.ADNC.ls, function(x) {
    x <- x %>%
        dplyr::filter(
            !is.na(LBD),
            !is.na(CVD),
            !is.na(LATE)
        )
    x
})
lapply(df.adj.ADNC.ls, dim)
DEA.adj3.ADNC.LMH.limma.res.numeric.ls <-
    diff_analysis(
        c("ADNC_LMH", "LBD", "CVD", "LATE", cov_names),
        df.adj.ADNC.ls,
        method = "limma",
        var_type = "numeric",
        NUM_COV = NUM_COV
    )

qs_savem(
    DEA.adj3.ADNC.LMH.limma.res.numeric.ls,
    file = paste0(RES_SUB_DIR, "DEA.adj.ADNC.LMH.limma.res.ls.QS2")
)

#----------------------------------------------------------------------------------------------------------------#
#-------------------------------------------------DONE-----------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------------#
#---------------------------------------------check results------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
if (TRUE) {
    qs_readm(paste0(RES_SUB_DIR, "DEA.limma.res.ls.QS2"))
    qs_readm(paste0(RES_SUB_DIR, "DEA.adj.ADNC.LMH.limma.res.ls.QS2"))

    All.sites <- mapply(
        function(a, b, c, d, e) {
            dfm <- topTable(m, coef = "A_beta_0_3", number = Inf)
            colnames(dfm) <- paste("a.score", colnames(dfm), sep = ".")
            dfn <- topTable(n, coef = "B_NFT_braak_0_3", number = Inf)
            colnames(dfn) <- paste("b.score", colnames(dfn), sep = ".")
            dfg <- topTable(g, coef = "ADNC_LMH", number = Inf)
            colnames(dfg) <- paste("ADNC.LMH.num", colnames(dfg), sep = ".")
            dfi <- topTable(i, coef = "Braak_NFT_stage", number = Inf)
            colnames(dfi) <- paste("braak.num", colnames(dfi), sep = ".")

            dfl <- topTable(l, coef = "ADNC_LMH", number = Inf)
            colnames(dfl) <- paste(
                "adj3.ADNC.LMH.num",
                colnames(dfl),
                sep = "."
            )

            if (
                sd(c(nrow(dfa), nrow(dfb), nrow(dfc), nrow(dfd), nrow(dfe))) !=
                    0
            ) {
                msg <- sprintf(
                    "The number of rows in the data frames is not the same: %d %d %d %d %d ",
                    nrow(dfa),
                    nrow(dfb),
                    nrow(dfc),
                    nrow(dfd),
                    nrow(dfe)
                )
                message(msg)
            }

            df.ls <- list(dfa, dfb, dfc, dfd, dfe)
            merged_df <- Reduce(
                function(x, y) {
                    x <- merge(x, y, by = "row.names", all = TRUE)
                    rownames(x) <- x$Row.names
                    x <- x[, -1]
                    x
                },
                df.ls
            )
            merged_df
        },

        DEA.a.score.limma.res.numeric.ls,
        DEA.b.score.limma.res.numeric.ls,
        DEA.ADNC.LMH.limma.res.numeric.ls,
        DEA.braak.limma.res.numeric.ls,
        DEA.adj3.ADNC.LMH.limma.res.numeric.ls,
        SIMPLIFY = F
    )
    lapply(All.sites, dim)
    param_combinations <- expand.grid(
        y = c(
            "a.score",
            "b.score",
            "ADNC.LMH.num",
            "braak.num",
            "adj3.ADNC.LMH.num"
        ),
        x = seq_along(All.sites)
    )

    # 使用 mapply 处理组合后的参数
    nothing.ls <- mapply(
        function(x, y) {
            p.col.name <- paste(y, "adj.P.Val", sep = ".")
            x.name <- names(All.sites)[x]

            df <- All.sites[[x]]
            idx <- which(colnames(df) == p.col.name)
            df <- df[df[, idx] < 0.05, ]
            df <- df[order(df[, idx]), ]

            str.print <- sprintf(
                "%s, %s < 0.05, %d sites",
                x.name,
                p.col.name,
                nrow(df)
            )
            message(str.print)
        },
        param_combinations$x,
        param_combinations$y,
        SIMPLIFY = F
    )
    rm(nothing.ls)
    gc()

    All.sites <- lapply(
        All.sites,
        get_gene_name,
        PTM = "phospho",
        just_gene = F,
        "Modified sequence",
        "Protein description",
        "Localization probability"
    )
    All.sites <- lapply(All.sites, function(x) {
        x$loc_on_MAPT_8 <- sapply(
            x$for_volcano,
            from_longest_isof_to_clinical_isof
        )
        x
    })
    qs_savem(All.sites, file = paste0(RES_SUB_DIR, "All.sites.limma.QS2"))
    qs_readm(paste0(RES_SUB_DIR, "All.sites.limma.QS2"))
    All.sites <- lapply(All.sites, function(x) {
        x <- x %>% dplyr::select(-ends_with("t"), -ends_with("B"))
        x <- x %>%
            dplyr::select(
                loc_on_MAPT_8,
                for_volcano,
                `Gene name`,
                `Modified sequence`,
                `Protein description`,
                starts_with("ADNC.LMH.num"),
                starts_with("braak.num"),
                everything()
            ) %>%
            dplyr::arrange(ADNC.LMH.num.adj.P.Val)
        x
    })

    head(All.sites[[1]])
    lapply(All.sites, ncol) # 55
    All.sites[[1]] %>% dplyr::filter(ADNC.LMH.num.adj.P.Val < 0.05) %>% nrow()
    write.xlsx(
        All.sites[[1]],
        file = paste0(RES_SUB_DIR, "All.sites.All.traits.intensity.limma.xlsx"),
        rowNames = TRUE
    )
}
#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------DONE-------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
