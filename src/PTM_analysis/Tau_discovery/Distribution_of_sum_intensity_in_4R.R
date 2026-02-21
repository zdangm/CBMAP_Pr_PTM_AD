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
    library(preprocessCore)
    library(impute)
    library(RColorBrewer)
    library(ggrepel)
    library(purrr)
    library(ggridges)
    library(conflicted)
})
rm(list = ls())
conflict_prefer_all("dplyr")
source("helper_func.R")
RES_SUB_DIR <- "../Results/Distribution_of_sum_of_intensity_in_tau/"
log2.raw <- list(
    phospho = "phos/",
    ubiq = "ubiq/",
    ace = "ace/"
)
NUM_COV <- 56


i <- 1
for (i in seq_along(log2.raw)) {
    ptm <- names(log2.raw)[i]
    if (ptm == "ace") {
        qs_readm(paste0(
            log2.raw[[i]],
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
        ))
    } else {
        qs_readm(paste0(
            log2.raw[[i]],
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.QS2"
        ))
    }
    df.MTBR <- log1p_raw.all.info.ls[[1]] %>% select(contains("MAPT"))
    tmp.cn <- sapply(colnames(df.MTBR), from_longest_isof_to_clinical_isof)
    useful <- which((tmp.cn != "not found on MAPT-8") & !is.na(tmp.cn))
    tmp.cn[useful] <- str_split_i(tmp.cn[useful], pattern = "_", i = 2)
    tmp.cn[useful] <- cut(
        as.numeric(tmp.cn[useful]),
        breaks = c(0, 150, 243, 368, 441),
        labels = c("N-terminal", "PRR", "MTBR", "C-terminal"),
        right = TRUE
    )
    domain.ls <- c("N-terminal", "PRR", "MTBR", "C-terminal")
    # tmp.cn[useful] <- cut(
    #     as.numeric(tmp.cn[useful]),
    #     breaks = c(0, 44, 73, 102, 150, 197, 243, 368, 441),
    #     labels = c("N-terminal", "1N", "2N", " ", "Proline-rich-1", "Proline-rich-2", "4R", "C-terminal"),
    #     right = TRUE
    # )
    covs <- lapply(log1p_raw.all.info.ls, function(x) {
        x <- x[, -seq.int(1, ncol(x) - NUM_COV)]
        x
    })
    df.MTBR.bak <- df.MTBR
    p.ls <- list()
    j <- 1
    for (j in seq.int(4)) {
        mydomain <- as.character(j)

        df.MTBR <- df.MTBR.bak[, tmp.cn == mydomain, drop = FALSE]
        if (is.null(df.MTBR) | ncol(df.MTBR) == 0) {
            message(mydomain)
            next
        }
        df.MTBR <- cbind(df.MTBR, covs[[1]])
        df.MTBR.hc <- df.MTBR %>%
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
        df.MTBR.case <- df.MTBR %>%
            filter(
                ADNC == 1,
                OTHER == 0,
                diag_ALS == 0,
                diag_SCZ == 0,
                diag_epilepsy == 0,
                diag_others == 0
            )
        df.MTBR.ADNC <- rbind(df.MTBR.hc, df.MTBR.case)
        df.MTBR.sum.ADNC <- df.MTBR.ADNC %>%
            mutate(totals = rowSums(select(., contains("MAPT")))) %>%
            select(totals, everything())
        df.MTBR.sum.ADNC$ADNC_LMH <- factor(
            df.MTBR.sum.ADNC$ADNC_LMH,
            levels = c(0, 1, 2, 3),
            label = c("HC", "ADNC_L", "ADNC_M", "ADNC_H"),
            ordered = T
        )
        mycolors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(4)

        p <- ggplot(
            df.MTBR.sum.ADNC,
            aes(x = totals, y = ADNC_LMH, fill = ADNC_LMH)
        ) +
            geom_density_ridges(alpha = 0.4, size = 0.4) + #
            scale_x_continuous(expand = c(.05, .05)) +
            # scale_y_discrete(expand =  expand_scale(mult = c(0.05, .65))) +  # for supplementary figure
            scale_y_discrete(expand = expand_scale(mult = c(0.05, .5))) + # for main figure
            labs(
                title = sprintf("%s intensity of %s", ptm, domain.ls[j]),
                x = "Summation",
                y = "ADNC_LMH"
            ) +
            # labs(title = "", x = "Summation", y = "ADNC_LMH") +
            scale_fill_manual(
                values = c(
                    "ADNC_H" = mycolors[4],
                    "ADNC_M" = mycolors[3],
                    "ADNC_L" = mycolors[2],
                    "HC" = mycolors[1]
                )
            ) +
            theme_bw() + #
            theme(
                axis.title = element_blank(),
                legend.title = element_blank(),
                axis.text.y = element_blank(),
                panel.grid = element_blank(),
                title = element_text(size = 14, face = "bold"),
                legend.position = "bottom"
            )
        p.ls[[mydomain]] <- p
    }
    mygrob <- marrangeGrob(
        p.ls,
        ncol = 1,
        nrow = 1,
        top = NULL
    )
    pdf(sprintf("%s%s.pdf", RES_SUB_DIR, ptm), width = 5, height = 12.) # for main Figure
    # pdf(sprintf("%s%s.pdf", RES_SUB_DIR, ptm), width = 6, height = 6.) # for Supplementary Figure
    print(mygrob)
    dev.off()
}
