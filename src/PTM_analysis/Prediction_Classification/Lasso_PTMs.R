#----------------------------------------------------------------------------------------------------------------#
#-------------------------------------------2025-05-15, LASSO code-----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

suppressPackageStartupMessages({
    library(openxlsx)
    library(magrittr)
    library(caret)
    library(patchwork)
    library(limma)
    library(UpSetR)
    library(gridExtra)
    library(qs2)
    library(optparse)
    library(ggrepel)
    library(impute)
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(performance)
    library(purrr)
    library(data.table)
    library(tidyverse)
    library(pROC)
    library(glmnet)
})
rm(list = ls())


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

log2_raw.729.dirs.ls <- list(
    phospho = "phos/",
    ubiq = "ubiq/",
    ace = "ace/",
    rna = "rna/",
    protein = "pr/",
    all = NULL
)
params <- expand.grid(ptm = names(log2_raw.729.dirs.ls), no_L = c(TRUE, FALSE))
cat("ptm =", opt$params, "\n")
mypar <- opt$params
ptm <- params[mypar, "ptm"]
no_L <- params[mypar, "no_L"]
conflicted::conflict_prefer_all("dplyr")
RES_SUB_DIR <- "../Results/Lasso_PTMs/"
source(
    "helper_func.R"
)
NUM_COV <- 56


need_adj <- F
ptm.ls <- list(
    phospho = list(),
    ubiq = list(),
    ace = list(),
    rna = list(),
    protein = list(),
    all = list()
)
six.roc.ls <- list(
    phospho = list(),
    ubiq = list(),
    ace = list(),
    rna = list(),
    protein = list(),
    all = list()
)
six.fit.ls <- list(
    phospho = list(),
    ubiq = list(),
    ace = list(),
    rna = list(),
    protein = list(),
    all = list()
)
six.best.lam.ls <- list(
    phospho = list(),
    ubiq = list(),
    ace = list(),
    rna = list(),
    protein = list(),
    all = list()
)
# ii <- 1
# for (ii in seq_along(ptm.ls)) {
# ptm <- names(ptm.ls)[ii]
if (ptm == "ace") {
    qs_readm(paste0(
        log2_raw.729.dirs.ls[[3]],
        "mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
    ))
} else if (ptm %in% c("ubiq", "phospho")) {
    qs_readm(paste0(
        log2_raw.729.dirs.ls[[ptm]],
        "mis_log2_impt_sc_rint.729.all.info.ls.QS2"
    ))
} else if (ptm %in% c("rna", "protein")) {
    qs_readm(paste0(
        log2_raw.729.dirs.ls[[ptm]],
        "log2_raw.729.all.info.ls.QS2"
    ))
} else {
    qs_readm(paste0(
        log2_raw.729.dirs.ls[["ace"]],
        "mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
    ))
    ace <- log1p_raw.all.info.ls
    sample_ord <- ace$intensity$id
    qs_readm(paste0(
        log2_raw.729.dirs.ls[["rna"]],
        "log2_raw.729.all.info.ls.QS2"
    ))
    rna <- log1p_raw.all.info.ls
    sample_ord <- intersect(sample_ord, rna$intensity$id)
    qs_readm(paste0(
        log2_raw.729.dirs.ls[["rna"]],
        "log2_raw.729.all.info.ls.QS2"
    ))
    rna <- log1p_raw.all.info.ls
    rna <- lapply(rna, function(x) {
        x <- x[match(sample_ord, x$id), ]
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
        colnames(x) <- paste("rna", colnames(x), sep = ".")
        x
    })

    qs_readm(paste0(
        log2_raw.729.dirs.ls[["protein"]],
        "log2_raw.729.all.info.ls.QS2"
    ))
    protein <- log1p_raw.all.info.ls
    protein <- lapply(protein, function(x) {
        x <- x[match(sample_ord, x$id), ]
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
        colnames(x) <- paste("protein", colnames(x), sep = ".")
        x
    })

    ace <- lapply(ace, function(x) {
        x <- x[match(sample_ord, x$id), ]
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
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
        x <- x[, seq.int(1, ncol(x) - NUM_COV)]
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
        colnames(x)[seq.int(1, ncol(x) - NUM_COV)] <- paste(
            "ubiq",
            colnames(x)[seq.int(1, ncol(x) - NUM_COV)],
            sep = "."
        )
        x
    })
    set.seed(2025)
    log1p_raw.all.info.ls <- mapply(
        function(r, pr, a, p, u) {
            x <- cbind(r, pr, a, p, u)
            if (T) {
                random_cols <- sample(
                    colnames(x)[seq.int(1, ncol(x) - NUM_COV)],
                    size = 17698
                )
                random_cols <- c(
                    random_cols,
                    colnames(x)[-seq.int(1, ncol(x) - NUM_COV)]
                )
                x <- x %>% select(all_of(random_cols))
            }
            x
        },
        rna,
        protein,
        ace,
        phospho,
        ubiq,
        SIMPLIFY = F
    )
}
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


if (no_L) {
    df.ADNC.ls <- lapply(df.ADNC.ls, function(x) {
        x <- x %>%
            filter(
                ADNC_LMH != 1
            )
        x
    })
}
lapply(df.ADNC.ls, dim)
table(df.ADNC.ls[[1]]$ADNC)
df.ADNC.resids.ls <- df.ADNC.ls
if (need_adj) {
    covs_regout.ls <- diff_analysis(
        c("sex_male", "PMD", "RIN", "age"),
        df.ADNC.ls,
        method = "glm",
        var_type = "factor_unordered",
        NUM_COV = NUM_COV
    )

    if (FALSE) {
        resids.normality.p.ls <- list(intensity = NULL, no_intensity = NULL)
        for (i in seq_along(covs_regout.ls)) {
            lm.ls <- covs_regout.ls[[i]]
            p.vals <- sapply(lm.ls, function(alm) {
                shapiro.test(residuals(alm))$p.value
            })
            resids.normality.p.ls[[i]] <- p.vals
        }
        hist(log(resids.normality.p.ls$intensity))
        abline(v = log(0.05), col = "red", lwd = 2)
        plot(dis.p(df.ADNC.ls[[1]]))
        plot(covs_regout.ls[[1]][[4]])
    }

    for (i in seq_along(covs_regout.ls)) {
        lm.ls <- covs_regout.ls[[i]]
        df <- df.ADNC.resids.ls[[i]]
        cnt <- 0
        for (j in seq_along(lm.ls)) {
            probe_id <- names(lm.ls)[j]
            myres <- residuals(lm.ls[[j]])
            df <- cbind(df, myres)
            rnd <- j %/% 1000
            if (cnt == rnd) {
                next
            }
            str <- sprintf(
                "cbind resids: %d / %d in %s",
                rnd * 1000,
                length(lm.ls),
                names(covs_regout.ls)[i]
            )
            message(str)
            cnt <- rnd
        }
        colnames(df) <- names(lm.ls)
        df.ADNC.resids.ls[[i]] <- df
    }
    lapply(df.ADNC.resids.ls, dim)
    df.ADNC.ls <- df.ADNC.resids.ls
}


coefs.ls <- list(intensity = list(), no_intensity = list())
roc.ls <- list(intensity = list(), no_intensity = list())
fit.ls <- list(intensity = list(), no_intensity = list())
bestlam.ls <- list(intensity = list(), no_intensity = list())


for (j in seq.int(1, 10)) {
    myseed <- 2025 + j
    set.seed(myseed)
    name_id <- paste0("seed_", myseed)

    train_index <- createDataPartition(
        df.ADNC.ls[[1]]$ADNC_LMH,
        p = 0.7,
        list = FALSE
    )
    for (i in seq_along(df.ADNC.resids.ls)) {
        df.ADNC.resids <- df.ADNC.resids.ls[[i]]
        df.ADNC <- df.ADNC.ls[[i]]

        df.ADNC.train <- df.ADNC[train_index, ]
        df.ADNC.test <- df.ADNC[-train_index, ]
        # df.ADNC.resids.train <- cbind(df.ADNC.resids[train_index, ], age = df.ADNC.train$age, sex_male = df.ADNC.train$sex_male)
        # df.ADNC.resids.test <- cbind(df.ADNC.resids[-train_index, ], age = df.ADNC.test$age, sex_male = df.ADNC.test$sex_male)
        cv_lasso <- cv.glmnet(
            as.matrix(df.ADNC.train[, seq.int(
                1,
                ncol(df.ADNC.train) - NUM_COV
            )]),
            df.ADNC.train$ADNC,
            alpha = 1,
            family = "binomial"
        )

        fit.ls[[i]][[name_id]] <- cv_lasso
        bestlam.ls[[i]][[name_id]] <- cv_lasso$lambda.min
        message(sprintf(
            "SEED: %d, ptm %s %s  best lambda: %f",
            myseed,
            ptm,
            names(df.ADNC.resids.ls)[i],
            bestlam.ls[[i]][[name_id]]
        ))
        coefs.ls[[i]][[name_id]] <- as.data.frame(as.matrix(coef(cv_lasso))) %>%
            filter(s1 != 0)
        coefs.ls[[i]][[name_id]] <- coefs.ls[[i]][[name_id]] %>%
            rownames_to_column("Biomarker")
        predP <- predict(
            cv_lasso,
            newx = as.matrix(df.ADNC.test[, seq.int(
                1,
                ncol(df.ADNC.train) - NUM_COV
            )]),
            s = "lambda.min",
            type = "response"
        )
        roc_obj <- roc(df.ADNC.test$ADNC, as.vector(predP))
        roc.ls[[i]][[name_id]] <- roc_obj
    }
}
six.roc.ls[[ptm]] <- roc.ls[[1]]
six.fit.ls[[ptm]] <- fit.ls[[1]]
six.best.lam.ls[[ptm]] <- bestlam.ls[[1]]
# }
qs_savem(
    roc.ls,
    fit.ls,
    file = sprintf(
        "%s%s_%s.QS",
        RES_SUB_DIR,
        ptm,
        if_else(no_L, "no_L", "with_L")
    )
)


if (F) {
    # plot
    myfit.tb.all <- tibble(
        sensitivity = numeric(),
        specificity = numeric(),
        auc = numeric(),
        with_L = logical(),
        ptm = character(),
        roc_id = character()
    )
    i <- "protein"
    jj <- 1
    for (i in c("rna", "protein", "phospho", "ubiq", "ace", "all")) {
        qs_readm(sprintf("%s%s_%s.QS", RES_SUB_DIR, i, "with_L"))
        for (jj in seq_along(roc.ls[[1]])) {
            myauc <- auc(roc.ls[[1]][[jj]])
            myroc <- pROC::smooth(roc.ls[[1]][[jj]], method = "density")
            tmpdf <- data.frame(
                sensitivity = myroc$sensitivities,
                specificity = myroc$specificities
            )
            tmpdf$auc <- myauc
            tmpdf$with_L <- T
            tmpdf$ptm <- i
            tmpdf$roc_id <- names(roc.ls[[1]])[jj]
            myfit.tb.all <- rbind(myfit.tb.all, tmpdf)
        }
        qs_readm(sprintf("%s%s_%s.QS", RES_SUB_DIR, i, "no_L"))
        for (jj in seq_along(roc.ls[[1]])) {
            myauc <- auc(roc.ls[[1]][[jj]])
            myroc <- pROC::smooth(roc.ls[[1]][[jj]], method = "density")
            tmpdf <- data.frame(
                sensitivity = myroc$sensitivities,
                specificity = myroc$specificities
            )
            tmpdf$auc <- myauc
            tmpdf$with_L <- F
            tmpdf$ptm <- i
            tmpdf$roc_id <- names(roc.ls[[1]])[jj]
            myfit.tb.all <- rbind(myfit.tb.all, tmpdf)
        }
    }
    myfit.tb.all %>% head()

    dim(myfit.tb.all)
    auc_labels <- myfit.tb.all %>%
        filter(ptm != "all") %>%
        group_by(ptm, with_L, roc_id) %>%
        summarise(auc = unique(auc), .groups = "drop") %>%
        group_by(ptm, with_L) %>%
        summarise(mean_auc = mean(auc), .groups = "drop") %>%
        mutate(
            label = sprintf("Avg AUC = %.3f", mean_auc),
            x = 0.95,
            y = 0.05
        ) #
    auc_labels
    auc_labels$ptm <- factor(
        auc_labels$ptm,
        levels = c("rna", "protein", "phospho", "ubiq", "ace"),
        labels = c("RNA", "Protein", "Phospho", "Ubiq", "Ace")
    )
    myfit.tb.all$ptm <- factor(
        myfit.tb.all$ptm,
        levels = c("rna", "protein", "phospho", "ubiq", "ace"),
        labels = c("RNA", "Protein", "Phospho", "Ubiq", "Ace")
    )

    p <- myfit.tb.all %>%
        ggplot(aes(
            x = 1 - specificity,
            y = sensitivity,
            color = as.factor(roc_id)
        )) +
        geom_line(alpha = 0.6, linewidth = 0.8) +
        facet_grid(with_L ~ ptm) +
        geom_label(
            data = auc_labels,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = 1,
            vjust = 0,
            size = 6,
            fill = "white",
            color = "black",
            label.size = 0.4
        ) +
        geom_segment(
            aes(x = 0, y = 0, xend = 1, yend = 1),
            linetype = "dashed",
            color = "grey50"
        ) +
        labs(
            title = "",
            x = "1 - Specificity",
            y = "Sensitivity",
            color = "roc_id"
        ) +
        theme_bw() +
        theme(
            legend.position = "none",
            strip.text = element_text(size = 14)
        )
    p
    ggsave(
        sprintf("%sroc_grid_plot.pdf", RES_SUB_DIR),
        plot = p,
        width = 16,
        height = 7
    )
}
