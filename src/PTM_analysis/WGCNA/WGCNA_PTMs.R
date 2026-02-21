#----------------------------------------------------------------------------------------------------------------#
#-------------------------------------------2025-04-05, WGCNA code-----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggdendro)
    library(data.table)
    library(openxlsx)
    library(magrittr)
    library(patchwork)
    library(limma)
    library(UpSetR)
    library(gridExtra)
    library(qs2)
    library(optparse)
    library(impute)
    library(gtools)
    library(purrr)
    library(WGCNA)
    library(ggrepel)
    library(KEGG.db)
    library(clusterProfiler)
    library(cowplot)
})
rm(list = ls())

option_list <- list(
    make_option(
        c("--ptm"),
        type = "double",
        default = 3,
        help = "ptm",
        metavar = "number"
    ),
    make_option(
        c("--datpre"),
        type = "double",
        default = 1,
        help = "datpre",
        metavar = "number"
    ),
    make_option(
        c("--site2gene_method"),
        type = "double",
        default = 2,
        help = "site2gene_method",
        metavar = "number"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
cat("ptm =", opt$ptm, "\n")
cat("datpre =", opt$datpre, "\n")
cat("site2gene_method =", opt$site2gene_method, "\n")
ptm <- opt$ptm
datpre <- opt$datpre
site2gene_method <- opt$site2gene_method


conflict_prefer_all("dplyr")

source("helper_func.R")
NUM_COV <- 56
log2_raw.729.dirs.ls <- list(
    phospho = "../phos/",
    ubiq = "../ubiq/",
    ace = "../ace/"
)
data_preprocesses <-
    c(
        "log_impute_scale_INT",
        "log_impute_scale"
    )

site2gene_methods <-
    c(
        "max",
        "median",
        "most_variable",

        "mean",
        "sum"
    )

if (F) {
    ptm <- 2
    datpre <- 1
    site2gene_method <- 2
}

use_sig_sites <- FALSE
do_not_redistribution <- TRUE

RES_SUB_DIR <- "../Results/WGCNA/"
RES_SUB_DIR <- paste0(RES_SUB_DIR, names(log2_raw.729.dirs.ls)[ptm], "/")
RES_SUB_DIR <- paste0(RES_SUB_DIR, data_preprocesses[datpre], "/")
RES_SUB_DIR <- paste0(RES_SUB_DIR, site2gene_methods[site2gene_method], "/")

if (use_sig_sites) {
    RES_SUB_DIR <- paste0(RES_SUB_DIR, "sig_sites/")
}

if (do_not_redistribution) {
    RES_SUB_DIR <- paste0(RES_SUB_DIR, "do_not_redistribution/")
}

dir.exists(RES_SUB_DIR)


if (T) {
    pipeline_steps <- strsplit(data_preprocesses[datpre], "_")[[1]]
    # qs_readm(paste0(log2_raw.729.dirs.ls[[ptm]], "log2_raw.729.all.info.ls.QS2"))
    if (ptm == 3) {
        qs_readm(paste0(
            log2_raw.729.dirs.ls[[ptm]],
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
        ))
    } else {
        qs_readm(paste0(
            log2_raw.729.dirs.ls[[ptm]],
            "log_impute_scale_INT/mis_log2_impt_sc_rint.729.all.info.ls.QS2"
        ))
    }

    log1p_raw.all.info.ls <- lapply(log1p_raw.all.info.ls, function(x) {
        x$sex_male <- factor(x$sex_male, ordered = FALSE)
        x
    })
    lapply(log1p_raw.all.info.ls, dim)

    if (use_sig_sites) {
        qs_readm(paste0(
            log2_raw.729.dirs.ls[[ptm]],
            data_preprocesses[datpre],
            "/All.sites.limma.QS2"
        ))
        All.sites <- lapply(All.sites, function(x) {
            x <- x %>% filter(ADNC.LMH.num.P.Value < 0.05)
            x
        })
        for (i in seq_along(All.sites)) {
            sig.sites <- rownames(All.sites[[i]])
            dat <- log1p_raw.all.info.ls[[i]]
            cov_names <- colnames(dat)[seq.int(
                ncol(dat) - NUM_COV + 1,
                ncol(dat)
            )]
            log1p_raw.all.info.ls[[i]] <- dat %>%
                dplyr::select(any_of(sig.sites), any_of(cov_names))
        }
    }

    sample_info <- log1p_raw.all.info.ls[[1]][, seq.int(
        ncol(log1p_raw.all.info.ls[[1]]) - NUM_COV + 1,
        ncol(log1p_raw.all.info.ls[[1]])
    )]
    sample_info$sex_male <- as.numeric(sample_info$sex_male)
    colnames(sample_info)
    head(sample_info[, 1:10])
    df.braak.ls <- lapply(log1p_raw.all.info.ls, function(x) {
        df.braak <- x %>%
            filter(
                OTHER == 0,
                diag_ALS == 0,
                diag_SCZ == 0,
                diag_epilepsy == 0,
                diag_others == 0,
                !is.na(Braak_NFT_stage)
            )
        df.braak
    })

    to_print <- sprintf(
        "braak.intensity: %d * %d, braak.no_intensity: %d * %d",
        dim(df.braak.ls$intensity)[1],
        dim(df.braak.ls$intensity)[2],
        dim(df.braak.ls$no_intensity)[1],
        dim(df.braak.ls$no_intensity)[2]
    )
    message(to_print)
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
    to_print <- sprintf(
        "adnc_hc.intensity: %d * %d, adnc_hc.no_intensity: %d * %d",
        dim(df.hc.ls$intensity)[1],
        dim(df.hc.ls$intensity)[2],
        dim(df.hc.ls$no_intensity)[1],
        dim(df.hc.ls$no_intensity)[2]
    )
    message(to_print)
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
    to_print <- sprintf(
        "adnc_case.intensity: %d * %d, adnc_case.no_intensity: %d * %d",
        dim(df.ADNC.case.ls$intensity)[1],
        dim(df.ADNC.case.ls$intensity)[2],
        dim(df.ADNC.case.ls$no_intensity)[1],
        dim(df.ADNC.case.ls$no_intensity)[2]
    )
    message(to_print)
    df.ADNC.ls <- mapply(
        function(x, y) {
            rbind(x, y)
        },
        df.hc.ls,
        df.ADNC.case.ls,
        SIMPLIFY = F
    )
    to_print <- sprintf(
        "adnc.intensity: %d * %d, adnc.no_intensity: %d * %d",
        dim(df.ADNC.ls$intensity)[1],
        dim(df.ADNC.ls$intensity)[2],
        dim(df.ADNC.ls$no_intensity)[1],
        dim(df.ADNC.ls$no_intensity)[2]
    )
    message(to_print)

    rm(df.hc.ls, df.ADNC.case.ls, log1p_raw.all.info.ls)
    names(df.ADNC.ls) <- paste("lmh", names(df.ADNC.ls), sep = ".")
    names(df.braak.ls) <- paste("braak", names(df.braak.ls), sep = ".")
    df0.lmh.braak.ls <- list(
        lmh.intensity = df.ADNC.ls$lmh.intensity,
        lmh.no_intensity = df.ADNC.ls$lmh.no_intensity,
        braak.intensity = df.braak.ls$braak.intensity,
        braak.no_intensity = df.braak.ls$braak.no_intensity
    )
    rm(df.ADNC.ls, df.braak.ls)

    i <- 1
    df0.lmh.braak.longer.ls <- list()
    for (i in seq_along(df0.lmh.braak.ls)) {
        x <- df0.lmh.braak.ls[[i]]
        probes <- x[, seq.int(1, ncol(x) - NUM_COV)]
        covs <- x[, seq.int(ncol(x) - NUM_COV + 1, ncol(x))]
        probes <- probes %>% rownames_to_column("sample")
        setDT(probes)
        probes_long <- melt(
            probes,
            id.vars = "sample",
            variable.name = "var",
            value.name = "intensity"
        )[,
            c("pro", "pos", "AA") := tstrsplit(var, "_", fixed = TRUE)
        ][,
            site := paste(pos, AA, sep = "_")
        ][,
            .(sample, pro, site, intensity)
        ]
        probes_long <- probes_long %>% as.data.frame()
        df0.lmh.braak.longer.ls[[names(df0.lmh.braak.ls)[i]]] <- probes_long
    }

    lapply(df0.lmh.braak.longer.ls, dim)

    head(df0.lmh.braak.longer.ls[[1]])
    myreshape <- function(x, site2gene_method = "max") {
        setDT(x)
        if (site2gene_method == "max") {
            x <- x[,
                .(gene_intensity = max(intensity, na.rm = TRUE)),
                by = .(sample, pro)
            ]
        } else if (site2gene_method == "mean") {
            x <- x[,
                .(gene_intensity = mean(intensity, na.rm = TRUE)),
                by = .(sample, pro)
            ]
        } else if (site2gene_method == "median") {
            x <- x[,
                .(gene_intensity = median(intensity, na.rm = TRUE)),
                by = .(sample, pro)
            ]
        } else if (site2gene_method == "most_variable") {
            x <- x[,
                sd_intensity := sd(intensity, na.rm = TRUE),
                by = .(pro, site)
            ]
            x <- x[, .SD[which.max(sd_intensity)], by = .(pro, sample)]
            x[, sd_intensity := NULL]
            x[, gene_intensity := intensity]
        } else if (site2gene_method == "sum") {
            x <- x[,
                .(gene_intensity = sum(intensity, na.rm = TRUE)),
                by = .(sample, pro)
            ]
        }
        x_wide <- dcast(x, sample ~ pro, value.var = "gene_intensity")
        setDF(x_wide)
        x_wide <- x_wide %>% column_to_rownames("sample")
        x_wide
    }
    i <- 1
    for (i in seq_along(df0.lmh.braak.ls)) {
        df0.lmh.braak.ls[[i]] <- myreshape(
            df0.lmh.braak.longer.ls[[i]],
            site2gene_methods[site2gene_method]
        )
    }
    for (i in seq_along(df0.lmh.braak.ls)) {
        df <- df0.lmh.braak.ls[[i]]
        df <- df %>% rownames_to_column("id")
        df0.lmh.braak.ls[[i]] <- df
    }

    lapply(df0.lmh.braak.ls, dim)
    # write_csv(df0.lmh.braak.ls$braak.intensity, "../tmp_for_brain_bank/phos_braak_mdian.csv")
    # write_csv(df0.lmh.braak.ls$braak.intensity, "../tmp_for_brain_bank/ubiq_braak_mdian.csv")
    # write_csv(df0.lmh.braak.ls$braak.intensity, "../tmp_for_brain_bank/ace_braak_mdian.csv")
    #----------------------------------------------------------------------------------------------------------------#
    #-------------------------------------------------begin----------------------------------------------------------#
    #----------------------------------------------------------------------------------------------------------------#
    mygsg <- function(x) {
        gsg <- goodSamplesGenes(x, verbose = 3)
        if (!gsg$allOK) {
            # If not all genes are good, remove the offending genes from the data.
            if (sum(!gsg$goodGenes) > 0) {
                message(paste(
                    "Removing genes:",
                    paste(names(x)[!gsg$goodGenes], collapse = ", ")
                ))
            }
            if (sum(!gsg$goodSamples) > 0) {
                message(paste(
                    "Removing samples:",
                    paste(rownames(x)[!gsg$goodSamples], collapse = ", ")
                ))
            }
            x <- x[gsg$goodSamples, gsg$goodGenes]
        }
        x
    }
    df0.lmh.braak.ls.bak <- df0.lmh.braak.ls
    for (i in seq_along(df0.lmh.braak.ls)) {
        df0.lmh.braak.ls[[i]] <- mygsg(df0.lmh.braak.ls[[i]])
    }
    lapply(df0.lmh.braak.ls, dim)

    obsvd_thds <- array(
        dim = c(3, 1, 5, 4),
        dimnames = list(
            Ptms = c("phospho", "ubiq", "ace"),
            Data_preprocesses = c(
                "log_impute_scale_INT"
            ),
            Site2gene_methods = c(
                "max",
                "sum",
                "mean",
                "median",
                "most_variable"
            ),
            Df0s = c(
                "lmh.intensity",
                "lmh.no_intensity",
                "braak.intensity",
                "braak.no_intensity"
            )
        )
    )

    obsvd_thds["phospho", "log_impute_scale_INT", "most_variable", ] <- c(
        110,
        54,
        110,
        50
    )
    obsvd_thds["phospho", "log_impute_scale_INT", "mean", ] <- c(70, 36, 70, 35)
    obsvd_thds["phospho", "log_impute_scale_INT", "max", ] <- c(80, 45, 76, 45)
    obsvd_thds["phospho", "log_impute_scale_INT", "median", ] <- c(
        79,
        35,
        85,
        34.5
    )
    obsvd_thds["phospho", "log_impute_scale_INT", "sum", ] <- c(
        375,
        180,
        400,
        180
    )

    obsvd_thds["ubiq", "log_impute_scale_INT", "most_variable", ] <- c(
        110,
        55,
        110,
        51
    )
    obsvd_thds["ubiq", "log_impute_scale_INT", "mean", ] <- c(75, 45, 80, 43)
    obsvd_thds["ubiq", "log_impute_scale_INT", "max", ] <- c(85, 50, 90, 50)
    obsvd_thds["ubiq", "log_impute_scale_INT", "median", ] <- c(80, 45, 90, 40)
    obsvd_thds["ubiq", "log_impute_scale_INT", "sum", ] <- c(500, 250, 500, 270)

    obsvd_thds["ace", "log_impute_scale_INT", "most_variable", ] <- c(
        63,
        30,
        63,
        25
    )
    obsvd_thds["ace", "log_impute_scale_INT", "mean", ] <- c(45, 25, 50, 25)
    obsvd_thds["ace", "log_impute_scale_INT", "max", ] <- c(50, 23, 52, 23)
    obsvd_thds["ace", "log_impute_scale_INT", "median", ] <- c(50, 20, 45, 20)
    obsvd_thds["ace", "log_impute_scale_INT", "sum", ] <- c(160, 65, 145, 65)

    sts <- list()
    pdf(file = paste0(RES_SUB_DIR, "sampleTree.pdf"), width = 20, height = 16)
    for (i in seq_along(df0.lmh.braak.ls)) {
        sampleTree <- hclust(dist(df0.lmh.braak.ls[[i]]), method = "average")
        sts[[names(df0.lmh.braak.ls)[i]]] <- sampleTree
        plot(
            sampleTree,
            main = paste0(
                "Sample clustering to detect outliers ",
                names(df0.lmh.braak.ls)[i]
            ),
            sub = "",
            xlab = "",
            cex = 0.6
        )
        if (!is.null(obsvd_thds)) {
            abline(
                h = obsvd_thds[
                    ptm,
                    datpre,
                    site2gene_method,
                    names(df0.lmh.braak.ls)[i]
                ],
                col = "red"
            )
        }
    }
    dev.off()
    qs_savem(sts, obsvd_thds, file = paste0(RES_SUB_DIR, "sampleTree.ls.QS2"))

    qs_readm(paste0(RES_SUB_DIR, "sampleTree.ls.QS2"))
    # df.lmh.braak.ls <- list()
    # for (i in seq_along(sts)) {
    #     sampleTree <- sts[[i]]
    #     clust <- cutreeStatic(
    #         sampleTree,
    #         cutHeight = obsvd_thds[ptm, datpre, site2gene_method, names(sts)[i]],
    #         minSize = 10
    #         )
    #     table(clust)
    #     # clust 1 contains the samples we want to keep.
    #     keepSamples <- (clust == 1)
    #     df.lmh.braak.ls[[names(sts)[i]]] <- df0.lmh.braak.ls[[i]][keepSamples, ]
    # }
    df.lmh.braak.ls <- df0.lmh.braak.ls
    datTraits.ls <- list()
    for (i in seq_along(df.lmh.braak.ls)) {
        datTraits.ls[[names(df.lmh.braak.ls)[i]]] <-
            sample_info[
                match(rownames(df.lmh.braak.ls[[i]]), sample_info$id),
            ] %>%
            dplyr::select(where(is.numeric))
        message(dim(datTraits.ls[[i]]))
    }

    pdf(
        file = paste0(RES_SUB_DIR, "Sample_dendrogram_and_trait_heatmap.pdf"),
        width = 20,
        height = 16
    )
    for (i in seq_along(datTraits.ls)) {
        sampleTree2 <- hclust(dist(df.lmh.braak.ls[[i]]), method = "average")
        datTraits <- datTraits.ls[[i]]
        traitColors <- numbers2colors(datTraits, signed = FALSE)
        plotDendroAndColors(
            sampleTree2,
            traitColors,
            groupLabels = names(datTraits),
            main = paste0(
                "Sample dendrogram and trait heatmap",
                names(df.lmh.braak.ls)[i]
            )
        )
    }
    dev.off()

    powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
    sft.ls <- list()
    power.ls <- list()
    enableWGCNAThreads(3)
    i <- 1
    pdf(
        file = paste0(RES_SUB_DIR, "Scale_free_topology_fit_index.pdf"),
        width = 10,
        height = 5
    )
    par(mfrow = c(1, 2))
    for (i in seq_along(df.lmh.braak.ls)) {
        datExpr <- df.lmh.braak.ls[[i]]
        sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        sft.ls[[names(df.lmh.braak.ls)[i]]] <- sft
        sft <- sft.ls[[names(df.lmh.braak.ls)[i]]]
        sft_data <- data.frame(
            Power = sft$fitIndices[, 1],
            R2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
            mean_connectivity = sft$fitIndices[, 5]
        )
        power_selected <- sft_data$Power[which(
            sft_data$R2 > 0.8 & sft_data$mean_connectivity < 100
        )[1]]

        str2print <- sprintf(
            "The selected power for %s is %d",
            names(df.lmh.braak.ls)[i],
            power_selected
        )
        message(str2print)
        # power_selected <- power.ls[[names(df.lmh.braak.ls)[i]]] # to delete
        power.ls[[names(df.lmh.braak.ls)[i]]] <- power_selected
        # Scale-free topology fit index as a function of the soft-thresholding power
        plot(
            sft$fitIndices[, 1],
            -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
            xlab = "Soft Threshold (power)",
            ylab = "Scale Free Topology Model Fit,signed R^2",
            type = "n",
            ylim = c(0, 1),
            # main = paste0("Scale independence")  # to delete
            main = paste0("Scale independence ", names(df.lmh.braak.ls)[i])
        )
        text(
            sft$fitIndices[, 1],
            -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
            labels = powers,
            cex = 0.9,
            col = "red"
        )
        abline(
            h = sft_data$R2[power_selected],
            v = power_selected,
            col = "red",
            lty = "dashed"
        )
        abline(h = 0.8, col = "red")
        # this line corresponds to using an R^2 cut-off of h
        # Mean connectivity as a function of the soft-thresholding power
        plot(
            sft$fitIndices[, 1],
            sft$fitIndices[, 5],
            xlab = "Soft Threshold (power)",
            ylab = "Mean Connectivity",
            type = "n",
            # main = paste0("Mean connectivity") # to delete
            main = paste0("Mean connectivity ", names(df.lmh.braak.ls)[i])
        )
        text(
            sft$fitIndices[, 1],
            sft$fitIndices[, 5],
            labels = powers,
            cex = 0.9,
            col = "red"
        )
        abline(
            h = sft_data$mean_connectivity[power_selected],
            v = power_selected,
            col = "red",
            lty = "dashed"
        )
        abline(h = 100, col = "red")
    }
    dev.off()
    qs_savem(
        sft.ls,
        power.ls,
        file = paste0(RES_SUB_DIR, "sft_and_power.ls.QS2")
    )

    df.lmh.braak.ls$lmh.no_intensity <- NULL
    df.lmh.braak.ls$braak.no_intensity <- NULL
    datTraits.ls$lmh.no_intensity <- NULL
    datTraits.ls$braak.no_intensity <- NULL
    power.ls$lmh.no_intensity <- NULL
    power.ls$braak.no_intensity <- NULL
    qs_savem(
        df.lmh.braak.ls,
        datTraits.ls,
        file = paste0(RES_SUB_DIR, "datExpr_and_datTraits.QS2")
    )

    qs_readm(paste0(RES_SUB_DIR, "sft_and_power.ls.QS2"))
    qs_readm(paste0(RES_SUB_DIR, "datExpr_and_datTraits.QS2"))
}


if (T) {
    qs_readm(paste0(RES_SUB_DIR, "sft_and_power.ls.QS2"))
    qs_readm(paste0(RES_SUB_DIR, "datExpr_and_datTraits.QS2"))
    lapply(df.lmh.braak.ls, dim)
    net.ls <- list()
    for (i in seq_along(df.lmh.braak.ls)) {
        datExpr <- df.lmh.braak.ls[[i]]
        net <- blockwiseModules(
            datExpr,
            power = power.ls[[names(df.lmh.braak.ls)[i]]],
            corType = "bicor",
            networkType = "signed",
            TOMType = "signed",
            deepSplit = 4,
            pamRespectsDendro = F,
            pamStage = F,
            minModuleSize = 30,
            mergeCutHeight = 0.07,
            numericLabels = TRUE,
            verbose = 3,

            maxBlockSize = 5000,
            reassignThreshold = 0,
            saveTOMs = TRUE,
            saveTOMFileBase = paste0(
                RES_SUB_DIR,
                names(df.lmh.braak.ls)[i],
                "-TOM"
            ),
            randomSeed = 2025,
            nThreads = 1 #
        )
        net.ls[[names(df.lmh.braak.ls)[i]]] <- net
    }

    qs_savem(net.ls, file = paste0(RES_SUB_DIR, "net.ls.QS2"))
}

if (T) {
    MEs.ls <- list()
    kME.ls <- list()
    moduleColors.ls <- list()
    color_table.ls <- list()
    hub_gene.ls <- list()
    qs_readm(paste0(RES_SUB_DIR, "net.ls.QS2"))
    qs_readm(paste0(RES_SUB_DIR, "datExpr_and_datTraits.QS2"))
    qs_readm(paste0(RES_SUB_DIR, "sft_and_power.ls.QS2"))

    i <- 1

    for (i in seq_along(net.ls)) {
        print(i)
        data <- df.lmh.braak.ls[[names(net.ls)[i]]]
        net <- net.ls[[i]]

        MEs <- moduleEigengenes(data, net$colors)$eigengenes
        kME <- signedKME(data, MEs, corFnc = "bicor")
        newColors <- net$colors

        for (module in colnames(kME)) {
            moduleNumber <- as.numeric(gsub("kME", "", module))
            inModule <- (newColors == moduleNumber)
            newColors[inModule & (kME[, module] < 0.3)] <- 0
        }
        if (do_not_redistribution == F) {
            maxIterations <- 2
            while (T) {
                changed <- 0
                grey_genes <- which(newColors == 0)
                if (length(grey_genes) == 0) {
                    break
                }

                for (ii in grey_genes) {
                    max_kME <- max(kME[ii, ], na.rm = TRUE)
                    if (max_kME > 0.3) {
                        bestModule <- which.max(kME[ii, ])
                        bestModuleNumber <- as.numeric(gsub(
                            "kME",
                            "",
                            colnames(kME)[bestModule]
                        ))
                        newColors[ii] <- bestModuleNumber
                        changed <- changed + (bestModuleNumber != 0)
                    }
                }
                if (changed == 0) break
            }
            net$colors <- newColors
            net.ls[[i]]$colors <- newColors
            MEs <- moduleEigengenes(data, net$colors)$eigengenes
            kME <- signedKME(data, MEs, corFnc = "bicor")
        }

        table(net$colors)

        moduleColors <- labels2colors(net$colors)
        table(moduleColors)
        Module_dendrogram <- sprintf(
            "%s%s_%s.pdf",
            RES_SUB_DIR,
            names(net.ls)[i],
            "Module_dendrogram"
        )
        pdf(file = Module_dendrogram, width = 6, height = 4)
        for (ii in seq_along(net$dendrograms)) {
            plotDendroAndColors(
                net$dendrograms[[ii]],
                moduleColors[net$blockGenes[[ii]]],
                "Modulecolors",
                dendroLabels = FALSE,
                hang = 0.03,
                addGuide = TRUE,
                guideHang = 0.05
            )
        }
        dev.off()

        color_table <- data.frame(
            ModuleNumber = sort(unique(net$colors)),
            ModuleColor = labels2colors(sort(unique(net$colors)))
        )

        kWithin <- intramodularConnectivity(
            adjacency(data, power = power.ls[[names(net.ls)[i]]]),
            net$colors
        )

        hub_pro <- chooseTopHubInEachModule(
            data,
            colorh = moduleColors,
            power = power.ls[[names(net.ls)[i]]],
            type = "signed"
        )
        # hub_gene <- pro_50_raw$gene_name[match(hub_pro, pro_50_raw$protein)]
        # hub_gene_modules <- moduleColors[hub_pro]

        hub_gene_module_info <- data.frame(
            Hub_Gene_Index = hub_pro,
            Module_Color = names(hub_pro),
            row.names = hub_pro
        )
        # hub_gene_module_info <- get_gene_name(hub_gene_module_info, PTM = names(log2_raw.729.dirs.ls)[ptm], just_gene = T)
        # hub_gene_module_info <- hub_gene_module_info %>% rename(Hub_Gene = `Gene name`)
        hub_gene_module_info$Hub_Gene <- str_split_i(
            hub_gene_module_info$Hub_Gene_Index,
            "\\.",
            1
        )
        MEs.ls[[names(net.ls)[i]]] <- MEs
        kME.ls[[names(net.ls)[i]]] <- kME
        moduleColors.ls[[names(net.ls)[i]]] <- moduleColors
        color_table.ls[[names(net.ls)[i]]] <- color_table
        hub_gene.ls[[names(net.ls)[i]]] <- hub_gene_module_info
    }

    legend_p_file_name <- sprintf("%s%s.pdf", RES_SUB_DIR, "legend_p")
    pdf(legend_p_file_name)
    for (i in seq_along(color_table.ls)) {
        color_table <- color_table.ls[[i]]

        max_per_col <- 8
        n <- nrow(color_table)

        if (n <= max_per_col) {
            color_table <- color_table %>%
                mutate(
                    y = max_per_col - row_number() + 1 #
                )

            p <- ggplot(color_table, aes(x = 1, y = y)) +
                geom_point(
                    aes(fill = ModuleColor),
                    shape = 21,
                    size = 10,
                    color = "black"
                ) +
                geom_text(
                    aes(x = 1.2, label = ModuleNumber),
                    hjust = 0,
                    size = 6
                ) +
                scale_fill_identity() +
                coord_cartesian(xlim = c(0.9, 2)) +
                theme_void() +
                theme(plot.margin = margin(10, 20, 10, 10))
        } else {
            color_table <- color_table %>%
                arrange(ModuleNumber) %>%
                mutate(
                    index = row_number(),
                    col = (index - 1) %/% max_per_col + 1,
                    y = max_per_col - ((index - 1) %% max_per_col)
                )

            color_table$ModuleNumber <- factor(
                color_table$ModuleNumber,
                levels = unique(color_table$ModuleNumber)
            )

            p <- ggplot(color_table, aes(x = 1, y = y)) +
                geom_point(
                    aes(fill = ModuleColor),
                    shape = 21,
                    size = 10,
                    color = "black"
                ) +
                geom_text(
                    aes(x = 1.2, label = ModuleNumber),
                    hjust = 0,
                    size = 6
                ) +
                scale_fill_identity() +
                facet_wrap(~col, nrow = 1, scales = "free_y") +
                coord_cartesian(xlim = c(0.9, 2)) +
                theme_void() +
                theme(
                    strip.text = element_blank(),
                    plot.margin = margin(10, 20, 10, 10)
                )
        }

        print(p)
    }
    dev.off()

    qs_savem(
        net.ls,
        MEs.ls,
        kME.ls,
        moduleColors.ls,
        color_table.ls,
        hub_gene.ls,
        file = paste0(
            RES_SUB_DIR,
            "MEs_kME_moduleColors_color_table_hub_gene.ls.QS2"
        )
    )
}

if (T) {
    qs_readm(paste0(
        RES_SUB_DIR,
        "MEs_kME_moduleColors_color_table_hub_gene.ls.QS2"
    ))

    pdf(
        file = paste0(RES_SUB_DIR, "hub_gene_barplot.pdf"),
        width = 6.4,
        height = 5
    )
    for (i in seq_along(net.ls)) {
        net <- net.ls[[i]]
        hub_gene <- hub_gene.ls[[names(net.ls)[i]]]
        mod_col <- data.frame(
            label = net$colors,
            color = labels2colors(net$colors),
            row.names = names(net$colors)
        )

        color_freq <- as.data.frame(table(mod_col$color))
        colnames(color_freq) <- c("Color", "Frequency")

        hub_genes_long <- data.frame(
            Color = hub_gene$Module_Color,
            Hub_Gene = hub_gene$Hub_Gene
        )

        color_freq <- color_freq %>%
            left_join(hub_genes_long, by = "Color") %>%
            arrange(desc(Frequency)) %>%
            mutate(
                Color = factor(Color, levels = unique(Color), ordered = TRUE)
            )
        matches <- data.frame(
            idx = paste0("ME", seq.int(nrow(color_freq)) - 1),
            Color = labels2colors(seq.int(nrow(color_freq)) - 1)
        )
        color_freq <- color_freq %>%
            left_join(matches, by = "Color")
        color_freq$Color <- factor(
            color_freq$Color,
            levels = matches$Color,
            ordered = TRUE
        )
        color_freq$idx <- factor(
            color_freq$idx,
            levels = matches$idx,
            ordered = TRUE
        )
        p <- ggplot(color_freq, aes(x = idx, y = Frequency)) +
            geom_bar(
                aes(fill = Color),
                stat = "identity",
                width = 0.6,
                color = "black",
                linewidth = 0.5
            ) +
            geom_text(
                aes(label = Frequency),
                hjust = -0.5,
                size = 5
                # family = "Times New Roman"
            ) +
            scale_fill_manual(
                values = as.character(color_freq$Color),
                labels = color_freq$Hub_Gene,
                name = "Hub"
            ) +
            scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
            coord_flip() +
            labs(
                x = "",
                y = "",
                title = "Module size"
                # title = paste0(names(net.ls)[i], ": Number of probes in each module with Hub Probe")
            ) +
            guides(fill = guide_legend(reverse = TRUE)) + #
            theme_bw() +
            theme(
                axis.text.x = element_text(size = 14, face = "bold"),
                axis.text.y = element_text(size = 14, face = "bold"),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.title = element_text(
                    size = 14,
                    face = "bold",
                    hjust = 0.5
                ),
                panel.grid = element_blank()
            )
        print(p)
    }
    dev.off()

    myenrich <- function(gene, module_color, n_probe, n_gene, idx) {
        GO_database <- "org.Hs.eg.db"
        KEGG_database <- "hsa"
        GO.BP <- enrichGO(
            gene$ENTREZID,
            OrgDb = GO_database,
            keyType = "ENTREZID",
            ont = "BP", # ontology = "ALL" means including Biological Process,Cellular Component,Mollecular Function
            readable = TRUE
        )
        title <- sprintf(
            "%s: %s probes of %s genes in module %s, %s",
            "BP",
            n_probe,
            n_gene,
            module_color,
            idx
        )

        GO.CC <- enrichGO(
            gene$ENTREZID,
            OrgDb = GO_database,
            keyType = "ENTREZID",
            ont = "CC", # ontology = "ALL" means including Biological Process,Cellular Component,Mollecular Function
            readable = TRUE
        )
        title <- sprintf(
            "%s: %s probes of %s genes in module %s, %s",
            "CC",
            n_probe,
            n_gene,
            module_color,
            idx
        )

        GO.MF <- enrichGO(
            gene$ENTREZID,
            OrgDb = GO_database,
            keyType = "ENTREZID",
            ont = "MF", # ontology = "ALL" means including Biological Process,Cellular Component,Mollecular Function
            readable = TRUE
        ) #
        title <- sprintf(
            "%s: %s probes of %s genes in module %s, %s",
            "MF",
            n_probe,
            n_gene,
            module_color,
            idx
        )

        tryCatch(
            {
                KEGG <- enrichKEGG(
                    gene$ENTREZID,
                    organism = KEGG_database,
                    use_internal_data = T ###
                )
            },
            error = function(e) {
                message("Error in enrichKEGG: ", e)
            }
        )
        title <- sprintf(
            "%s: %s probes of %s genes in module %s, %s",
            "KEGG",
            n_probe,
            n_gene,
            module_color,
            idx
        )

        # plot_BP <- dotplot(GO.BP, font.size = 10) + ggtitle(title)
        # plot_CC <- dotplot(GO.CC, font.size = 10) + ggtitle(title)
        # plot_MF <- dotplot(GO.MF, font.size = 10) + ggtitle(title)
        # plot_KEGG <- dotplot(KEGG) + ggtitle(title)

        # p <- list(BP = plot_BP, CC = plot_CC, MF = plot_MF, KEGG = plot_KEGG)
        df <- list(BP = GO.BP, CC = GO.CC, MF = GO.MF, KEGG = KEGG)
        return(list(p = NULL, df = df))
    }

    myenrich2 <- function(gene, idx) {
        GO_database <- "org.Hs.eg.db"
        KEGG_database <- "hsa"
        GO <- enrichGO(
            gene$ENTREZID,
            OrgDb = GO_database,
            keyType = "ENTREZID",
            ont = "ALL", # ontology = "ALL" means including Biological Process,Cellular Component,Mollecular Function
            readable = TRUE
        )
        tryCatch(
            {
                KEGG <- enrichKEGG(
                    gene$ENTREZID,
                    organism = KEGG_database,
                    use_internal_data = T ###
                )
            },
            error = function(e) {
                message("Error in enrichKEGG: ", e)
            }
        )
        GO.bak <- GO
        GO.p <- module_highlight(
            go_enrich = GO,
            GOtopN = 10
        )
        KEGG.p <- module_highlight(
            kegg_enrich = KEGG,
            KEGGtopN = 25
        )
        all.p <- module_highlight(
            go_enrich = GO,
            kegg_enrich = KEGG,
            GOtopN = 10,
            KEGGtopN = 10
        )

        return(list(GO = GO.p, KEGG = KEGG.p, all = all.p))
    }

    p.ls <- list()

    mod_col.ls <- list()
    for (i in seq_along(net.ls)) {
        net <- net.ls[[i]]
        mod_col <- data.frame(
            label = net$colors,
            color = labels2colors(net$colors),
            row.names = names(net$colors)
        )
        mod_col$`Gene name` <- str_split_i(rownames(mod_col), "\\.", 1)
        #
        # mod_col %>% rownames_to_column("id") %>% arrange(label) %>% fwrite(file = paste0(RES_SUB_DIR, "gene_and_module.txt"), sep = "\t")
        mod_col.ls[[names(net.ls)[i]]] <- mod_col
        for (j in unique(mod_col$label)) {
            # for (j in js) {
            message(j)
            if (j == 0) {
                next
            }
            label_col <- labels2colors(j)
            message(label_col)
            mod_col_sub <- mod_col %>% filter(label == j)
            gene_EID <- bitr(
                mod_col_sub$`Gene name`,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = "org.Hs.eg.db"
            )
            gene <- gene_EID %>%
                dplyr::select(ENTREZID) %>%
                unique()
            idx <- mod_col_sub$label[1]
            n_gene <- nrow(gene)

            keggo_sub_dir <- "keggo_plots_dfs/"
            dir.create(
                paste0(RES_SUB_DIR, keggo_sub_dir),
                showWarnings = T,
                recursive = TRUE
            )

            # result <- myenrich(gene, module_color = j, n_probe = nrow(mod_col_sub), idx = idx, n_gene = n_gene)

            result2 <- tryCatch(
                {
                    myenrich2(gene, idx = idx)
                },
                error = function(e) {
                    message("Error in myenrich2: ", e)
                    return(NULL)
                }
            )
            if (is.null(result2)) {
                message(
                    "Skipping module ",
                    labels2colors(j),
                    " in network ",
                    names(net.ls)[i],
                    " due to error in myenrich2."
                )
                next
            }
            p.ls[[j]] <- result2$all #
            p_file_name2 <- sprintf(
                "%s%s%s_%s_beautiful_keggo.pdf",
                RES_SUB_DIR,
                keggo_sub_dir,
                names(net.ls)[i],
                labels2colors(j)
            )
            pdf(p_file_name2, width = 7, height = 12)
            tryCatch(
                {
                    print(result2$all)
                },
                error = function(e) {
                    message("empty plot: ", e)
                }
            )
            dev.off()

            df.BP <- result$df$BP@result
            df.CC <- result$df$CC@result
            df.MF <- result$df$MF@result
            df.KEGG <- result$df$KEGG@result
            df.KEGG <- kegg_results_fixing(df.KEGG)

            wb <- createWorkbook()
            addWorksheet(wb, "GO_BP")
            writeData(wb, sheet = "GO_BP", df.BP)

            addWorksheet(wb, "GO_CC")
            writeData(wb, sheet = "GO_CC", df.CC)

            addWorksheet(wb, "GO_MF")
            writeData(wb, sheet = "GO_MF", df.MF)
            addWorksheet(wb, "KEGG")
            df.KEGG <- kegg_results_fixing(df.KEGG)
            writeData(wb, sheet = "KEGG", df.KEGG)
            df_file_name <- sprintf(
                "%s%s%s_%s.xlsx",
                RES_SUB_DIR,
                keggo_sub_dir,
                names(net.ls)[i],
                labels2colors(j)
            )
            saveWorkbook(wb, file = df_file_name, overwrite = TRUE)
        }
        # p.ls.sub <- list()
        # for (ii in seq_along(p.ls)) {
        #     if (!is.null(p.ls[[ii]])) {
        #         p.ls.sub[[as.character(ii)]] <- p.ls[[ii]]
        #     }
        # }
        # pdf("tmp.pdf", width = 20, height = 24)
        # p.ls.sub[[1]] + p.ls.sub[[2]] + p.ls.sub[[3]] + p.ls.sub[[4]] + p.ls.sub[[5]] + p.ls.sub[[6]] + p.ls.sub[[7]] + p.ls.sub[[8]] + plot_layout(ncol = 4, guides = "collect")
        # dev.off()
    }
}

if (T) {
    qs_readm(paste0(
        RES_SUB_DIR,
        "MEs_kME_moduleColors_color_table_hub_gene.ls.QS2"
    ))
    qs_readm(paste0(RES_SUB_DIR, "datExpr_and_datTraits.QS2"))
    i <- 1
    trait_table.ls <- list()
    for (i in seq_along(df.lmh.braak.ls)) {
        datExpr <- df.lmh.braak.ls[[i]]
        datTraits <- datTraits.ls[[names(df.lmh.braak.ls)[i]]]
        datTraits <-
            datTraits %>%
            dplyr::select(
                # -(year:RIN),
                # -OTHER,
                # -(diag_ALS:diag_others),
                # -hc,
                # -CVD,
                # -(ADNC_L:ADNC_H),
                # -(CVD_atherosclerosis:CVD_CAA),
                # -starts_with("ARTAG")
                ADNC_LMH,
                age,
                sex_male,
                PMD,
                RIN,
                Braak_NFT_stage,
                A_beta_0_3,
                B_NFT_braak_0_3,
                C_cerad_0_3,
                PART,
                LATE,
                ARTAG,
                LBD,
                CVD
            ) %>%
            rownames_to_column("id")
        nGenes <- ncol(datExpr)
        nSamples <- nrow(datExpr)
        moduleColors <- moduleColors.ls[[names(df.lmh.braak.ls)[i]]]
        MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
        MEs <- orderMEs(MEs0)
        MEs <- t(MEs)

        trait_adnc <- datTraits %>% filter(!is.na(ADNC_LMH))
        me_adnc <- MEs[, colnames(MEs) %in% trait_adnc$id]
        me_adnc <- me_adnc[, match(trait_adnc$id, colnames(me_adnc))]
        design <- model.matrix(
            ~ ADNC_LMH + age + sex_male + PMD + RIN,
            data = trait_adnc
        )
        fit <- lmFit(me_adnc, design)
        fit <- eBayes(fit)
        adnc_num_lim_re <- topTable(
            fit,
            coef = "ADNC_LMH",
            number = Inf,
            adjust.method = "BH"
        )
        adnc_num_lim_re$me <- rownames(adnc_num_lim_re)
        adnc_num_lim_re$trait <- "ADNC"
        # head(adnc_num_lim_re)

        age_lim_re <- topTable(
            fit,
            coef = "age",
            number = Inf,
            adjust.method = "BH"
        )
        age_lim_re$me <- rownames(age_lim_re)
        # head(age_lim_re)
        age_lim_re$trait <- 'Age'

        # trait_braak <- datTraits %>% filter(!is.na(Braak_NFT_stage))
        # dim(trait_braak)
        # me_braak <- MEs[, colnames(MEs) %in% trait_braak$id]
        # dim(me_braak)
        # me_braak <- me_braak[, match(trait_braak$id, colnames(me_braak))]
        # dim(me_braak)
        # design <- model.matrix(~ Braak_NFT_stage + age + sex_male + PMD + RIN, data = trait_braak)
        # fit <- lmFit(me_braak, design)
        # fit <- eBayes(fit)
        # braak_num_lim_re <- topTable(fit, coef = "Braak_NFT_stage", number = Inf, adjust.method = "BH")
        # braak_num_lim_re$me <- rownames(braak_num_lim_re)
        # dim(braak_num_lim_re)
        # braak_num_lim_re$trait <- "braak"

        trait_ascore <- datTraits %>% filter(!is.na(A_beta_0_3))
        # dim(trait_ascore)
        me_ascore <- MEs[, colnames(MEs) %in% trait_ascore$id]
        # dim(me_ascore)
        me_ascore <- me_ascore[, match(trait_ascore$id, colnames(me_ascore))]
        # dim(me_ascore)
        design <- model.matrix(
            ~ A_beta_0_3 + age + sex_male + PMD + RIN,
            data = trait_ascore
        )
        fit <- lmFit(me_ascore, design)
        fit <- eBayes(fit)
        ascore_num_lim_re <- topTable(
            fit,
            coef = "A_beta_0_3",
            number = Inf,
            adjust.method = "BH"
        )
        ascore_num_lim_re$me <- rownames(ascore_num_lim_re)
        # dim(ascore_num_lim_re)
        ascore_num_lim_re$trait <- "Abeta"

        trait_bscore <- datTraits %>% filter(!is.na(B_NFT_braak_0_3))
        # dim(trait_bscore)
        me_bscore <- MEs[, colnames(MEs) %in% trait_bscore$id]
        # dim(me_bscore)
        me_bscore <- me_bscore[, match(trait_bscore$id, colnames(me_bscore))]
        # dim(me_bscore)
        design <- model.matrix(
            ~ B_NFT_braak_0_3 + age + sex_male + PMD + RIN,
            data = trait_bscore
        )
        fit <- lmFit(me_bscore, design)
        fit <- eBayes(fit)
        bscore_num_lim_re <- topTable(
            fit,
            coef = "B_NFT_braak_0_3",
            number = Inf,
            adjust.method = "BH"
        )
        bscore_num_lim_re$me <- rownames(bscore_num_lim_re)
        # dim(bscore_num_lim_re)
        bscore_num_lim_re$trait <- "Braak"

        trait_cscore <- datTraits %>% filter(!is.na(C_cerad_0_3))
        # dim(trait_cscore)
        me_cscore <- MEs[, colnames(MEs) %in% trait_cscore$id]
        # dim(me_cscore)
        me_cscore <- me_cscore[, match(trait_cscore$id, colnames(me_cscore))]
        # dim(me_cscore)
        design <- model.matrix(
            ~ C_cerad_0_3 + age + sex_male + PMD + RIN,
            data = trait_cscore
        )
        fit <- lmFit(me_cscore, design)
        fit <- eBayes(fit)
        cscore_num_lim_re <- topTable(
            fit,
            coef = "C_cerad_0_3",
            number = Inf,
            adjust.method = "BH"
        )
        cscore_num_lim_re$me <- rownames(cscore_num_lim_re)
        dim(cscore_num_lim_re)
        cscore_num_lim_re$trait <- "Cerad"

        trait_part <- datTraits %>% filter(!is.na(PART))
        dim(trait_part)
        me_part <- MEs[, colnames(MEs) %in% trait_part$id]
        dim(me_part)
        me_part <- me_part[, match(trait_part$id, colnames(me_part))]
        dim(me_part)
        design <- model.matrix(
            ~ PART + age + sex_male + PMD + RIN,
            data = trait_part
        )
        fit <- lmFit(me_part, design)
        fit <- eBayes(fit)
        part_num_lim_re <- topTable(
            fit,
            coef = "PART",
            number = Inf,
            adjust.method = "BH"
        )
        part_num_lim_re$me <- rownames(part_num_lim_re)
        dim(part_num_lim_re)
        part_num_lim_re$trait <- "PART"

        trait_late <- datTraits %>% filter(!is.na(LATE))
        dim(trait_late)
        me_late <- MEs[, colnames(MEs) %in% trait_late$id]
        dim(me_late)
        me_late <- me_late[, match(trait_late$id, colnames(me_late))]
        dim(me_late)
        design <- model.matrix(
            ~ LATE + age + sex_male + PMD + RIN,
            data = trait_late
        )
        fit <- lmFit(me_late, design)
        fit <- eBayes(fit)
        late_num_lim_re <- topTable(
            fit,
            coef = "LATE",
            number = Inf,
            adjust.method = "BH"
        )
        late_num_lim_re$me <- rownames(late_num_lim_re)
        dim(late_num_lim_re)
        late_num_lim_re$trait <- "LATE"

        trait_artag <- datTraits %>% filter(!is.na(ARTAG))
        dim(trait_artag)
        me_artag <- MEs[, colnames(MEs) %in% trait_artag$id]
        dim(me_artag)
        me_artag <- me_artag[, match(trait_artag$id, colnames(me_artag))]
        dim(me_artag)
        design <- model.matrix(
            ~ ARTAG + age + sex_male + PMD + RIN,
            data = trait_artag
        )
        fit <- lmFit(me_artag, design)
        fit <- eBayes(fit)
        artag_num_lim_re <- topTable(
            fit,
            coef = "ARTAG",
            number = Inf,
            adjust.method = "BH"
        )
        artag_num_lim_re$me <- rownames(artag_num_lim_re)
        dim(artag_num_lim_re)
        artag_num_lim_re$trait <- "ARTAG"

        trait_lbd <- datTraits %>% filter(!is.na(LBD))
        dim(trait_lbd)
        me_lbd <- MEs[, colnames(MEs) %in% trait_lbd$id]
        dim(me_lbd)
        me_lbd <- me_lbd[, match(trait_lbd$id, colnames(me_lbd))]
        dim(me_lbd)
        design <- model.matrix(
            ~ LBD + age + sex_male + PMD + RIN,
            data = trait_lbd
        )
        fit <- lmFit(me_lbd, design)
        fit <- eBayes(fit)
        lbd_num_lim_re <- topTable(
            fit,
            coef = "LBD",
            number = Inf,
            adjust.method = "BH"
        )
        lbd_num_lim_re$me <- rownames(lbd_num_lim_re)
        dim(lbd_num_lim_re)
        lbd_num_lim_re$trait <- "LBD"

        trait_cvd <- datTraits %>% filter(!is.na(CVD))
        dim(trait_cvd)
        me_cvd <- MEs[, colnames(MEs) %in% trait_cvd$id]
        dim(me_cvd)
        me_cvd <- me_cvd[, match(trait_cvd$id, colnames(me_cvd))]
        dim(me_cvd)
        design <- model.matrix(
            ~ CVD + age + sex_male + PMD + RIN,
            data = trait_cvd
        )
        fit <- lmFit(me_cvd, design)
        fit <- eBayes(fit)
        cvd_num_lim_re <- topTable(
            fit,
            coef = "CVD",
            number = Inf,
            adjust.method = "BH"
        )
        cvd_num_lim_re$me <- rownames(cvd_num_lim_re)
        dim(cvd_num_lim_re)
        cvd_num_lim_re$trait <- "CVD"

        trait_table <- rbind(
            adnc_num_lim_re,
            age_lim_re,
            # braak_num_lim_re,
            bscore_num_lim_re,
            ascore_num_lim_re,
            cscore_num_lim_re,
            part_num_lim_re,
            late_num_lim_re,
            artag_num_lim_re,
            lbd_num_lim_re,
            cvd_num_lim_re
        )
        trait_table.ls[[names(df.lmh.braak.ls)[i]]] <- trait_table
    }
    trait_table.ls <- lapply(trait_table.ls, function(x) {
        x$pmarker <- ifelse(x$adj.P.Val < 0.05, "*", "")
        x$me <- factor(x$me)
        x$trait <- factor(x$trait)
        x
    })
    trait_table.ls <- lapply(trait_table.ls, function(x) {
        num_col <- length(unique(x$me))
        matches <- paste0("ME", labels2colors(seq.int(num_col) - 1))
        matches <- data.frame(
            idx = paste0("ME", seq.int(num_col) - 1),
            me = matches
        )
        x <- x %>% left_join(matches, by = "me")
        x
    })
    p.ls <- list()
    i <- 1

    write.xlsx(trait_table.ls[[1]], paste0(RES_SUB_DIR, "trait_ME_assoc.xlsx"))

    for (i in seq_along(trait_table.ls)) {
        trait_table.ls[[i]]$trait <- factor(
            trait_table.ls[[i]]$trait,
            levels = c(
                "ADNC",
                "Age",
                "Abeta",
                "Braak",
                "Cerad",
                "PART",
                "LATE",
                "ARTAG",
                "LBD",
                "CVD"
            )
        )
        df <- trait_table.ls[[i]]
        df$idx <- factor(df$idx, levels = mixedsort(unique(df$idx)))
        df$trait <- factor(df$trait, levels = unique(df$trait))
        p <- ggplot(df, aes(x = idx, y = trait, fill = -log10(adj.P.Val))) +
            geom_tile(color = "white") +
            geom_text(
                aes(label = pmarker),
                color = "black",
                size = 9,
                vjust = 0.5,
                hjust = 0.5
            ) +
            scale_fill_gradient(
                low = "white",
                high = "red",
                name = bquote("-log"[10] * "(P)")
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(
                    angle = 45,
                    hjust = 1,
                    size = 14,
                    color = "black",
                    face = "bold"
                ),
                axis.text.y = element_text(
                    size = 14,
                    face = "bold",
                    color = "black"
                ),
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 8),
                panel.grid = element_blank()
            ) +
            labs(x = NULL, y = NULL, title = NULL)
        p.ls[[names(trait_table.ls)[i]]] <- p
    }
    pdf(
        file = paste0(RES_SUB_DIR, "Module_trait_relationships.pdf"),
        width = 6,
        height = 4
    )
    print(p.ls[[1]])
    print(p.ls[[2]])
    dev.off()
}
