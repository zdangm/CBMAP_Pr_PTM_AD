#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#

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
})
rm(list = ls())
source("helper_func.R")
conflicted::conflict_prefer_all("dplyr")
RES_SUB_DIR.ls <- list(
    phospho = "phos/",
    ubiq = "ubiq/",
    ace = "ace/"
)
NUM_COV <- 56

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

        # x <- x %>% dplyr::select(matches("^MAPT\\.|^SYNPR\\.|^APP\\."))
        x <- x %>% dplyr::select(matches("^MAPT\\."))

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
identical(rownames(MAPT.ls.1[[1]]), batch_info$`分析名`)
# TRUE
dim(batch_info)
# 729 4
# sample_info <- fread("/data/shared_data/China_Brain_MultiOmics/sample_information/final/CBMAP_sample_info_1187_final_20241101.csv")
sample_info <- fread(
    "/data/shared_data/China_Brain_MultiOmics/sample_information/final/CBMAP_sample_info_1187_final_20250523.csv"
)

sample_info$id <- if_else(
    sample_info$bank == "zju",
    paste0("A", sample_info$id),
    sample_info$id
)
sample_info <- sample_info[match(rownames(MAPT.ls.1[[1]]), sample_info$id), ]
identical(rownames(MAPT.ls.1[[1]]), sample_info$id)
# TRUE
dim(sample_info)
"PTB087" %in% sample_info$id
# TRUE
"repeat2" %in% sample_info$id
# FALSE
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
df.ADNC.ls.bak <- df.ADNC.ls
#----------------------------------------------------------------------------------------------------------------#
#-------------------------------------above is the same with Circle_heatmap.R------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
df.ADNC.ls <- lapply(df.ADNC.ls, function(x) {
    x <- x %>%
        dplyr::select(
            matches("^phospho|^ubiq|^ace"),
            ADNC_LMH,
            age,
            sex_male,
            RIN,
            PMD
        )
    x.covs <- x %>% dplyr::select(-matches("^phospho|^ubiq|^ace"))
    x.probes <- x %>% dplyr::select(matches("^phospho|^ubiq|^ace"))
    cn <- colnames(x.probes)
    ptms <- str_split_i(cn, "\\.", 1)
    sites <- paste0(str_split_i(cn, "\\.", 2), ".", str_split_i(cn, "\\.", 3))
    new_cn <- sapply(sites, from_longest_isof_to_clinical_isof)
    new_cn <- if_else(is.na(new_cn), sites, new_cn)
    idx <- new_cn == "not found on MAPT-8"
    colnames(x.probes) <- paste0(ptms, ".", new_cn)
    x.probes <- x.probes[, !idx]
    x <- cbind(x.probes, x.covs)
    x
})
lapply(df.ADNC.ls, dim)
colnames(df.ADNC.ls$intensity)

#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------except 262-----------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
cov_names <- c("sex_male", "age", "PMD", "RIN")
top1_cov_names <- c(cov_names, "phospho.Tau441_262_S")
top2_cov_names <- c(cov_names, "phospho.Tau441_262_S", "phospho.Tau441_289_S")
top3_cov_names <- c(
    cov_names,
    "phospho.Tau441_262_S",
    "phospho.Tau441_289_S",
    "phospho.Tau441_212_T"
)
top4_cov_names <- c(
    cov_names,
    "phospho.Tau441_262_S",
    "phospho.Tau441_289_S",
    "phospho.Tau441_212_T",
    "phospho.Tau441_217_T"
)

top1_cov_names <- c(cov_names, "ubiq.Tau441_311_K")
top2_cov_names <- c(cov_names, "ubiq.Tau441_311_K", "ubiq.Tau441_321_K")
top3_cov_names <- c(
    cov_names,
    "ubiq.Tau441_311_K",
    "ubiq.Tau441_321_K",
    "ubiq.Tau441_254_K"
)
top4_cov_names <- c(
    cov_names,
    "ubiq.Tau441_311_K",
    "ubiq.Tau441_321_K",
    "ubiq.Tau441_254_K",
    "ace.Tau441_274_K"
)
test_cov_names <- c(
    cov_names,
    "ubiq.Tau441_311_K",
    "ubiq.Tau441_321_K",
    "ubiq.Tau441_254_K",
    "ace.Tau441_274_K",
    "ubiq.Tau441_281_K"
)


test_cov_names <- c(cov_names, "phospho.Tau441_217_T")

pure_assoc_p <- list(
    top1 = list(),
    top2 = list(),
    top3 = list(),
    top4 = list(),
    test = list()
)
# pure_assoc_p <- list(phospho = list(), ubiq = list(), ace = list())
res <- data.frame(
    ADNC_LMH = df.ADNC.ls$intensity$ADNC_LMH
)
cov_names.ls <- list(
    top1 = top1_cov_names,
    top2 = top2_cov_names,
    top3 = top3_cov_names,
    top4 = top4_cov_names,
    test = test_cov_names
)
for (i in seq_along(df.ADNC.ls[[1]])) {
    for (j in seq_along(cov_names.ls)) {
        var_name <- colnames(df.ADNC.ls[[1]])[i]

        if (var_name %in% cov_names.ls[[j]]) {
            next
        }
        if (var_name == "ADNC_LMH") {
            next
        }
        message("Processing ", var_name)
        my_formula <- sprintf(
            "%s ~ %s",
            var_name,
            paste(cov_names.ls[[j]], collapse = " + ")
        )
        fit <- lm(as.formula(my_formula), data = df.ADNC.ls$intensity)
        res[[var_name]] <- resid(fit)
        my_formula <- sprintf("%s ~ ADNC_LMH", var_name)
        pure_assoc <- lm(as.formula(my_formula), data = res)
        pure_assoc_p[[j]][[var_name]] <- summary(pure_assoc)$coefficients[2, 4]
        message("P-value for ", var_name, " is ", pure_assoc_p[[j]][[var_name]])
    }
}
length(pure_assoc_p)
pure_assoc_p <- lapply(pure_assoc_p, function(x) {
    name_x <- names(x)
    idx <- str_detect(name_x, ".*Tau441.*")
    x <- x[idx]
    x
})

pure_assoc_p[[1]] %>% unlist() %>% p.adjust(method = "BH") %>% sort()
pure_assoc_p[[1]] %>% unlist() %>% sort() %>% names()
pure_assoc_p[[1]] %>% unlist() %>% sort()
pure_assoc_p[[2]] %>% unlist() %>% p.adjust(method = "BH") %>% sort()
pure_assoc_p[[2]] %>% unlist() %>% sort() %>% names()
pure_assoc_p[[2]] %>% unlist() %>% sort()
pure_assoc_p[[3]] %>% unlist() %>% p.adjust(method = "BH") %>% sort()
pure_assoc_p[[3]] %>% unlist() %>% sort() %>% names()
pure_assoc_p[[3]] %>% unlist() %>% sort()
pure_assoc_p[[4]] %>% unlist() %>% p.adjust(method = "BH") %>% sort()
pure_assoc_p[[4]] %>% unlist() %>% sort() %>% names()
pure_assoc_p[[4]] %>% unlist() %>% sort()
pure_assoc_p[[5]] %>% unlist() %>% p.adjust(method = "BH") %>% sort()
pure_assoc_p[[5]] %>% unlist() %>% sort() %>% names()
pure_assoc_p[[5]] %>% unlist() %>% sort()


qs_savem(
    pure_assoc_p,
    file = "pure_assoc_p.QS2"
)

qs_savem(
    pure_assoc_p,
    file = "sequential_assoc_p.QS2"
)
