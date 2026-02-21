# install.packages("/home/lik/Scripts/likun/Random_tasks/Other_source_data/KEGG.db_1.0.tar.gz", repos = NULL, type = "source")
library(KEGG.db)
mylog2 <- function(x, NUM_COV) {
    message(sprintf("Log2 transformation, %d covs", NUM_COV))
    x.covs <- NULL
    if (NUM_COV > 0) {
        x.covs <- x[, seq.int(ncol(x) - NUM_COV + 1, ncol(x))]
    }
    x.probes <- x[, seq.int(1, ncol(x) - NUM_COV)]
    x.probes <- log2(x.probes + 1)
    if (!is.null(x.covs)) {
        x.probes <- cbind(x.probes, x.covs)
    }
    x.probes
}
myimputeknn <- function(x, NUM_COV) {
    message(sprintf("Impute missing values using kNN, %d covs", NUM_COV))
    x.covs <- NULL
    if (NUM_COV > 0) {
        x.covs <- x[, seq.int(ncol(x) - NUM_COV + 1, ncol(x))]
    }
    x.mat <- x[, seq.int(1, ncol(x) - NUM_COV)] %>% as.matrix()
    x.mat.t <- t(x.mat)
    # seems like impute.knn() requires the data should be probes x samples
    # so we need to transpose the data
    # x.mat.t.imptd <- impute.knn(x.mat.t, k = 10, rng.seed = 2025)$data
    x.mat.t.imptd <- suppressMessages({
        set.seed(2025)
        impute.knn(x.mat.t, k = 10, rng.seed = 2025)$data
    })
    x.mat.imptd <- t(x.mat.t.imptd)
    x.imptd <- x.mat.imptd %>% as.data.frame(check.names = F)
    if (!is.null(x.covs)) {
        x.imptd <- merge(x.imptd, x.covs, by = "row.names", all = FALSE)
        rownames(x.imptd) <- x.imptd$Row.names
        x.imptd <- x.imptd[, -1]
    }
    x.imptd
}
mymedian <- function(x, NUM_COV) {
    message(sprintf("Median centering, %d", NUM_COV))
    x[, seq.int(1, ncol(x) - NUM_COV)] <- apply(x[, seq.int(1, ncol(x) - NUM_COV)], 1, function(y) {
        y - median(y, na.rm = T)
    }) %>%
        t() %>% # deadly important transpose
        as.data.frame(check.names = F)
    x
}


myimputedmedian <- function(x, NUM_COV) {
    message(sprintf("imputed-Median centering, then fill back NA, %d", NUM_COV))
    x.covs <- NULL
    if (NUM_COV > 0) {
        x.covs <- x[, seq.int(ncol(x) - NUM_COV + 1, ncol(x))]
    }
    x.mat <- x[, seq.int(1, ncol(x) - NUM_COV)] %>% as.matrix()
    bool.x.mat <- is.na(x.mat)
    x.mat.t <- t(x.mat)
    # seems like impute.knn() requires the data should be probes x samples
    # so we need to transpose the data
    x.mat.t.imptd <- suppressMessages({
        set.seed(2025)
        impute.knn(x.mat.t, k = 10, rng.seed = 2025)$data
    })
    x.mat.imptd <- t(x.mat.t.imptd)
    x.mat.imptd <- apply(x.mat.imptd, 1, function(y) {
        y - median(y, na.rm = T)
    }) %>% t()
    x.mat.imptd[bool.x.mat] <- NA
    x.imptd <- x.mat.imptd %>% as.data.frame(check.names = F)
    if (!is.null(x.covs)) {
        x.imptd <- merge(x.imptd, x.covs, by = "row.names", all = FALSE)
        rownames(x.imptd) <- x.imptd$Row.names
        x.imptd <- x.imptd[, -1]
    }
    x.imptd
}

mynormalize <- function(x, NUM_COV = 56) {
    message(sprintf("mean 0, sd 1, %d covs", NUM_COV))
    x[, seq.int(1, ncol(x) - NUM_COV)] <- apply(x[, seq.int(1, ncol(x) - NUM_COV)], 1, function(y) {
        scale(y, na.rm = T)
    }) %>%
        t() %>% # deadly important transpose
        as.data.frame(check.names = F)
    x
}
myINT <- function(x, NUM_COV = 56) {
    message(sprintf("INT, need to be done after imputation! %d covs", NUM_COV))
    x[, seq.int(1, ncol(x) - NUM_COV)] <- apply(x[, seq.int(1, ncol(x) - NUM_COV)], 2, function(y) {
        qqnorm(y, plot.it = F)$x
    }) %>% as.data.frame(check.names = F)
    x
}

preprocess_in_DEA <- list(mylog2, myimputeknn, mymedian, myimputedmedian, mynormalize, myINT)
names(preprocess_in_DEA) <- c("log", "impute", "scale", "imscale", "normalize", "INT")




kegg_results_fixing <- function(df.kegg) {
    if (nrow(df.kegg) == 0) {
      return(df.kegg)
    }
    kegg_gene_ids <- df.kegg$geneID
    kegg_gene_ids <- strsplit(kegg_gene_ids, "/")
    # j <- 140
    for (j in seq_along(kegg_gene_ids)) {
        if (length(kegg_gene_ids[[j]]) == 1 & all(is.na(kegg_gene_ids[[j]]))) {
            kegg_gene_ids[[j]] <- "NA"
            next
        }
        suppressMessages({
            result <- bitr(
                kegg_gene_ids[[j]],
                fromType = "ENTREZID",
                toType = c("SYMBOL"),
                OrgDb = "org.Hs.eg.db",
                drop = T
            )
        })
        kegg_gene_ids[[j]] <- paste(result$SYMBOL, collapse = "/")
    }
    kegg_gene_ids <- data.frame(
        unlist(kegg_gene_ids),
        stringsAsFactors = FALSE,
        check.names = F
    )
    df.kegg$gene_name <- kegg_gene_ids[[1]]
    df.kegg
}












diff_analysis <- function(variables, data_list, method = c("glm", "limma"), var_type = c("factor_ordered", "factor_unordered", "numeric"), NUM_COV) {
    # 参数校验
    method <- match.arg(method)
    var_type <- match.arg(var_type)

    # 第一个变量是关注的变量
    target_var <- variables[1]

    # 构建设计矩阵
    build_design_matrix <- function(df, variables) {
        formula_str <- paste("~", paste(variables, collapse = " + ")) # 生成公式字符串
        model_matrix <- model.matrix(as.formula(formula_str), data = df)
        print(colnames(model_matrix))
        return(model_matrix)
    }

    # 差异分析函数
    perform_analysis <- function(df, method, var_type) {
        expression_matrix <- df[, -seq.int(ncol(df) - NUM_COV + 1, ncol(df)), drop = F] # 去除NUM_COV列
        covs <- df[, variables, drop = F]
        # 匹配目标变量的类型
        if (var_type == "factor_ordered") {
            covs[[target_var]] <- factor(covs[[target_var]], ordered = TRUE)
        } else if (var_type == "factor_unordered") {
            covs[[target_var]] <- factor(covs[[target_var]], ordered = FALSE)
        } else if (var_type == "numeric") {
            covs[[target_var]] <- as.numeric(covs[[target_var]])
        }
        if (method == "glm") {
            # 针对GLM的分析：目标变量为行、探针为列

            glm.fit.res <- list()
            progress <- 1
            for (i in seq_along(expression_matrix)) {
                fit <- glm(expression_matrix[, i] ~ ., data = covs, family = gaussian)
                phos_site <- colnames(expression_matrix)[i]
                glm.fit.res[[phos_site]] <- fit
                if (progress %% 100 == 0) {
                    message("Progress: ", progress, " / ", ncol(expression_matrix))
                }
                progress <- progress + 1
            }
            return(glm.fit.res)
        } else if (method == "limma") {
            # 针对limma的分析：样本为列、探针为行
            expression_matrix <- t(expression_matrix) # 转置
            design_matrix <- build_design_matrix(covs, variables)
            fit <- limma::lmFit(expression_matrix, design_matrix)
            fit <- limma::eBayes(fit)
            return(fit)
        }
    }

    # 遍历data_list并对每个数据框执行差异分析
    results <- lapply(data_list, function(df) {
        perform_analysis(df, method, var_type)
    })

    return(results)
}

draw_volcano_plot <- function(model.ls, var_name = NULL, fc_threshold = 0, p_threshold = 0.05, gene_name = NA, show_all = TRUE, PTM = NULL) {
    fc_threshold <- fc_threshold
    p_threshold <- p_threshold
    type <- NULL
    if (is.null(PTM)) {
        stop("please specify the PTM type")
    }
    if (!PTM %in% c("phospho", "ace", "ubiq")) {
        stop("PTM should be phospho, ace or ubiq")
    }

    if (all(sapply(model.ls, inherits, "glm"))) {
        df <- extract_glm_stats(model.ls, var_name)
        df <- get_gene_name(df, PTM = PTM)
        df <- df %>% rename(logFC = Estimate)
        type <- "glm"
    } else if (inherits(model.ls, "MArrayLM")) {
        df <- topTable(model.ls, coef = var_name, number = Inf)
        df <- get_gene_name(df, PTM = PTM)
        df$isofrom441 <- sapply(df$for_volcano, from_longest_isof_to_clinical_isof)
        df$isofrom441 <- if_else(!is.na(df$isofrom441), if_else(str_detect(df$isofrom441, "^Tau"), df$isofrom441, df$for_volcano), df$isofrom441)
        type <- "limma"
    } else {
        stop("model.ls should be a MArrayLM object or a list of glm objects")
    }
    df$logFC <- 2^(df$logFC) # exponen transformation
    # 注意，df里面是有NA的！！！也就导致了Direction中也有NA
    df$Direction <- ifelse(
        # abs(df$logFC) > fc_threshold & df$adj.P.Val < p_threshold,
        abs(log2(df$logFC)) > log2(fc_threshold) & df$P.Value < p_threshold,
        ifelse(df$logFC > 2^0, "Up", "Down"),
        "Not Significant"
    )
    df$neg_log10_p <- -log10(df$P.Value)

    # 选出需要标注的基因
    df$label_gene <- ""
    # 当 gene_name不为NA时，show_all参数才有意义。show_all==TRUE代表把gene_name里面的基因的所有sites都在火山图上标注
    num_sigs <- sum(df$adj.P.Val < p_threshold)
    num_up_sigs <- sum(df$adj.P.Val < p_threshold & df$logFC > 2^0)
    num_down_sigs <- sum(df$adj.P.Val < p_threshold & df$logFC < 2^0)
    if (all(!is.na(gene_name))) {
        if ("MAPT" %in% gene_name) {
            tau_rows <- grepl("^Tau441", df$isofrom441)
            sig_rows <- (df$Direction != "Not Significant") & (!is.na(df$Direction))

            df$label_gene[tau_rows & sig_rows] <- str_extract(df$isofrom441[tau_rows & sig_rows], "^[^_]+_[^_]+")
            if (show_all) {
                df$label_gene[tau_rows] <- str_extract(df$isofrom441[tau_rows], "^[^_]+_[^_]+")
            }
        }

        gene_name <- gene_name[gene_name != "MAPT"]
        if (length(gene_name) != 0) {
            gene_rows <- df$`Gene name` %in% gene_name

            sig_rows <- (df$Direction != "Not Significant") & (!is.na(df$Direction))

            df$label_gene[gene_rows & sig_rows] <- str_extract(df$for_volcano[gene_rows & sig_rows], "^[^_]+_[^_]+")
            if (show_all) {
                df$label_gene[gene_rows] <- str_extract(df$for_volcano[gene_rows], "^[^_]+_[^_]+")
            }
        }
    } else {
        sig_rows <- (df$Direction != "Not Significant") & (!is.na(df$Direction))
        df$label_gene[sig_rows] <- str_extract(df$for_volcano[sig_rows], "^[^_]+_[^_]+")
        if (any(str_detect(df$label_gene[sig_rows], "^MAPT"))) {
            tau_rows <- grepl("^Tau441", df$isofrom441)
            df$label_gene[tau_rows & sig_rows] <- str_extract(df$isofrom441[tau_rows & sig_rows], "^[^_]+_[^_]+")
        }
    }


    # 如果出现重复标注，只保留一个（防止重叠）
    df$label_gene <- ifelse(duplicated(df$label_gene), "", df$label_gene)

    # 开始绘图
    p <- ggplot(df, aes(x = logFC, y = neg_log10_p, color = Direction)) +
        geom_point(size = 2., alpha = 0.7) +
        scale_color_manual(values = c("Up" = "#D7263D", "Down" = "#1B9CFC", "Not Significant" = "gray80")) +
        geom_vline(xintercept = c(2^(-log2(fc_threshold)), fc_threshold), linetype = "dashed", color = "black") +
        geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "black") +
        geom_text_repel(
            data = df[df$label_gene != "", ],
            aes(label = label_gene),
            min.segment.length = 0,
            direction = "both", # 允许 x 和 y 方向移动
            force = 25, # 增加排斥力，让标签更分散
            force_pull = 0.1, # 减小吸引力，标签离点稍远
            point.padding = 1, # 标签与点的距离
            box.padding = 0.6, # 标签之间的距离
            max.overlaps = 15, # 允许更多重叠
            segment.size = 0.5, # 连接线粗细
            segment.alpha = 0.5, # 连接线透明度
            seed = 42, # 固定随机种子
            segment.color = "#676565", # 连接线颜色
            show.legend = FALSE
        ) +
        labs(
            x = "Fold Change",
            y = "-log10(P-value)",
            title = sprintf("Volcano Plot of %s\n%d sites fdr < 0.05 (%d +, %d -)", var_name, num_sigs, num_up_sigs, num_down_sigs),
            color = "Regulation"
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.position = "top"
        )

    return(list(df = df, p = p))
}



my_miami_plot <- function() {
    library(biomaRt)

    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    genes <- c("CAMK2G", "PSAP", "SYNPR", "TMEFF1", "TRIM47")

    annot <- getBM(
        attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
        filters = "hgnc_symbol",
        values = genes,
        mart = mart
    )
}


























from_longest_isof_to_clinical_isof <- function(from_loc) {
    # 将我们数据中的tau蛋白isoform的位置转换为临床中常用的tau蛋白isoform的位置
    if (is.na(from_loc)){
      return(from_loc)
    }
    from_loc <- as.character(from_loc)
    AA <- ""
    if (str_detect(from_loc[1], "_")) {
        arr <- str_split(from_loc, "_", simplify = TRUE)
        from_loc <- arr[1, 2]
        from_loc <- as.numeric(from_loc)
        Gene <- arr[1, 1]
        if (str_detect(Gene, "\\.")) {
            Gene <- str_split_i(Gene, "\\.", 1)
        }
        if (Gene != "MAPT") {
            return(NA)
        }
        AA <- arr[1, 3]
    } else {
        from_loc <- as.numeric(from_loc)
    }
    to_loc <- NA
    if (from_loc < 125) {
        to_loc <- from_loc
    } else if (from_loc > 375 & from_loc < 395) {
        to_loc <- from_loc - 251
    } else if (from_loc > 460) {
        to_loc <- from_loc - 317
    } else {
        to_loc <- "not found on MAPT-8"
        return(to_loc)
    }
    to_loc <- as.character(to_loc)
    if (AA != "") {
        to_loc <- paste0(to_loc, "_", AA)
    }
    to_loc <- paste0("Tau441_", to_loc)
    return(to_loc)
}




extract_glm_stats <- function(glm_list, var_name = NULL) {
    results <- t(sapply(glm_list, function(model) {
        coef_table <- summary(model)$coefficients
        if (!is.null(var_name)) {
            idx <- grep(var_name, rownames(coef_table))
            if (length(idx) == 0) stop("Variable not found in coefficients")
        } else {
            idx <- 2
        }
        c(
            t.value = coef_table[idx, "t value"],
            Estimate = coef_table[idx, "Estimate"],
            Std_err = coef_table[idx, "Std. Error"],
            Pr_t = coef_table[idx, "Pr(>|t|)"]
        )
    }))
    results <- data.frame(results, check.names = F)
    results$adj.P.Val <- p.adjust(results$Pr_t, method = "BH")
    results <- results %>% arrange(adj.P.Val)
    return(results)
}

get_gene_name <- function(df, PTM = NULL, just_gene = F, ...) {
    # df should be a data frame with phos-sites rownames: Protein accession_Position_Amino acid
    # ... is the extra columns you want to add from original intensity data frame
    intensity_ms_idnttfd_info <- NULL
    if (is.null(PTM)) {
        stop("please specify the PTM type")
    }
    if (PTM == "phospho") {
        # intensity_ms_idnttfd_info <- fread("/data/shared_data/China_Brain_MultiOmics/humanBrain_Phospho/XB07045B4DPST_post_mut_987samples_Dedupe_update_20250717/MS_identified_information.txt")
        intensity_ms_idnttfd_info <- fread("/data/shared_data/China_Brain_MultiOmics/humanBrain_Phospho/XB07045B4DPST_0311/XB07045B4DPST_mix+sample/L0G0/MS_identified_information.txt")
    } else if (PTM == "ubiq") {
        intensity_ms_idnttfd_info <- fread("/data/shared_data/China_Brain_MultiOmics/huamnBrain_ubiquitylation/XB07045B4DPUb_初分析报告释放/泛素化修饰组结果/MS_identified_information.txt")
    } else if (PTM == "ace") {
        intensity_ms_idnttfd_info <- fread("/data/shared_data/China_Brain_MultiOmics/huamnBrain_Acetylation/乙酰化修饰结果/MS_identified_information.txt")
    } else {
        stop("no such PTM type")
    }
    if (is.null(dim(intensity_ms_idnttfd_info))) {
        stop("intensity_ms_idnttfd_info is empty")
    }
    info <- intensity_ms_idnttfd_info[, c(1:7)]
    just_gene
    str_count(rownames(df)[1], "_")
    if (just_gene & (str_count(rownames(df)[1], "_") < 2)) {
        info <- info %>% dplyr::select(`Protein accession`, `Gene name`)
        info <- info[match(rownames(df), info$`Protein accession`), ]
        info <- info %>% dplyr::select(`Gene name`)
        df <- cbind(df, info)
        return(df)
    }
    info <- info %>% unite("id", "Protein accession", "Position", "Amino acid", sep = "_", remove = FALSE)
    info$id <- paste(info$`Gene name`, info$id, sep = ".")
    info <- info %>% unite("for_volcano", "Gene name", "Position", "Amino acid", sep = "_", remove = FALSE)
    fixed_cols <- c("id", "for_volcano", "Gene name")
    extra_cols <- c(...)
    selected_cols <- unique(c(fixed_cols, extra_cols))
    if (length(setdiff(selected_cols, colnames(info))) > 0) {
        stop("Some columns are not found in the info data frame.")
    }
    info <- info %>% dplyr::select(all_of(selected_cols))
    info <- info[match(rownames(df), info$id), ]
    df <- cbind(df, info)
    df
}



dis.p <- function(x, seed = 2025, NUM_COV) {
    # x should be sample * (phos_site + covariates)
    # 从数据框中随便抽几个位点出来看分布
    p.ls <- list()
    set.seed(seed)
    y <- sample(seq.int(1, ncol(x) - NUM_COV), 12)
    y <- colnames(x)[y]
    p.ls <- lapply(y, function(i) {
        # 提取每一列的数据
        data_col <- x[, i]

        # 进行正态性检验，获取p-value
        test_result <- shapiro.test(data_col)
        p_value <- test_result$p.value

        # 创建直方图并注释p-value
        p <- ggplot(data.frame(data_col, check.names = F), aes(x = data_col)) +
            geom_histogram(bins = 100) +
            ggtitle(i) +
            theme_minimal() +
            annotate("text",
                x = Inf, y = Inf, label = paste("p =", round(p_value, 4)),
                hjust = 1.1, vjust = 2, size = 5, color = "red"
            )

        p
    })
    tmp.p <- do.call(gridExtra::arrangeGrob, c(p.ls, list(ncol = 3, nrow = 4)))
    tmp.p
}

plot_phospho_distribution <- function(df, site_col, group_col) {
    if (!str_detect(site_col, "_")) {
        tmp <- str_split_i(site_col, "-", i = 1)
        tmp <- str_split(tmp, ";", simplify = TRUE)
        pos <- str_extract(tmp[1, 2], "\\d+")
        aa <- str_extract(tmp[1, 2], "[A-Z]")
        site_col <- sprintf("%s_%s_%s", tmp[1, 1], pos, aa)
    }
    # 提取分组列和位点列，移除缺失值
    df_subset <- df %>%
        dplyr::select(all_of(group_col), all_of(site_col)) %>%
        filter(!is.na(.data[[group_col]]), !is.na(.data[[site_col]])) %>%
        mutate(!!group_col := as.factor(.data[[group_col]])) # 转为因子
    # 核密度图
    p_density <- ggplot(df_subset, aes(x = .data[[site_col]], fill = .data[[group_col]])) +
        geom_density(alpha = 0.5) +
        labs(
            title = paste("Density Plot of", site_col, "by", group_col),
            x = "Signal Intensity",
            y = "Density"
        ) +
        theme_minimal() +
        theme(legend.position = "top")
    # 箱线图
    p_box <- ggplot(df_subset, aes(x = .data[[group_col]], y = .data[[site_col]], fill = .data[[group_col]])) +
        geom_boxplot() +
        labs(
            title = paste("Boxplot of", site_col, "by", group_col),
            x = group_col,
            y = "Signal Intensity"
        ) +
        theme_minimal() +
        theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1)) +
        stat_compare_means(method = "t.test", aes(group = .data[[group_col]]), label = "p.signif", label.y = max(df_subset[[site_col]], na.rm = TRUE) * 1.1)
    p.ls <- list(p_density, p_box)
    p.ls
}




enrichment_gene <- function(
    model.ls,
    var_name = NULL,
    PTM = NULL,
    enrich_method = c("GO", "KEGG"),
    fdr.cut.be4.keggo = 0.05,
    n_show_term = 25) {
    # 注意，相同的模型，limma里面的logFC和glm里面的Estimate是一样的意义
    enrich_method <- match.arg(enrich_method)
    if (is.null(PTM)) {
        stop("please specify the PTM type")
    }
    if (!PTM %in% c("phospho", "ace", "ubiq")) {
        stop("PTM should be phospho, ace or ubiq")
    }
    # Gene symbol to entrez id
    GO_database <- "org.Hs.eg.db"
    KEGG_database <- "hsa"
    if (inherits(model.ls, "MArrayLM")) {
        df <- topTable(model.ls, coef = var_name, number = Inf)
        df <- get_gene_name(df, PTM = PTM)
        type <- "limma"
        gene <- bitr(df$`Gene name`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = GO_database)
        df <- df %>% inner_join(gene, by = c("Gene name" = "SYMBOL"))
        df <- df %>% filter(adj.P.Val < fdr.cut.be4.keggo)
        gene <- df %>%
            dplyr::select(ENTREZID) %>%
            unique()
    } else if (is.character(model.ls)) {
        # 如果直接就是基因名 字符串向量
        gene <- bitr(model.ls, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = GO_database)
        gene <- unique(gene)
    } else {
        stop("model.ls should be a MArrayLM object or a list of glm objects or just a vector of gene names")
    }

    if (enrich_method == "GO") {
        GO <- enrichGO(
            gene$ENTREZID,
            OrgDb = GO_database,
            keyType = "ENTREZID",
            ont = "ALL", # ontology = "ALL" means including Biological Process,Cellular Component,Mollecular Function
            readable = TRUE
        ) # 设定是否将GO ID转换为Term
        df <- GO
    } else if (enrich_method == "KEGG") {
        KEGG <- enrichKEGG(
            gene$ENTREZID,
            organism = KEGG_database,
            use_internal_data = T ### 注意！！
        )
        df <- KEGG
    }
    return(df)
}
if (F) {
    GO.1 <- read.xlsx("/share/home/lik/Scripts/likun/Random_tasks/Results/WGCNA/ace/log_impute_scale_INT/median/do_not_redistribution/keggo_plots_dfs/lmh.intensity_pink.xlsx", sheet = 1)
    GO.2 <- read.xlsx("/share/home/lik/Scripts/likun/Random_tasks/Results/WGCNA/ace/log_impute_scale_INT/median/do_not_redistribution/keggo_plots_dfs/lmh.intensity_pink.xlsx", sheet = 2)
    GO.3 <- read.xlsx("/share/home/lik/Scripts/likun/Random_tasks/Results/WGCNA/ace/log_impute_scale_INT/median/do_not_redistribution/keggo_plots_dfs/lmh.intensity_pink.xlsx", sheet = 3)
    GO.1$ONTOLOGY <- "BP"
    GO.2$ONTOLOGY <- "CC"
    GO.3$ONTOLOGY <- "MF"
    GO <- rbind(GO.1, GO.2, GO.3)
    KEGG <- read.xlsx("/share/home/lik/Scripts/likun/Random_tasks/Results/WGCNA/ace/log_impute_scale_INT/median/do_not_redistribution/keggo_plots_dfs/lmh.intensity_pink.xlsx", sheet = 4)
    KEGG$ONTOLOGY <- "KEGG"
    GOtopN <- 15
    KEGGtopN <- 15

    GO <- GO %>%
        group_by(ONTOLOGY) %>%
        arrange(p.adjust, -Count) %>%
        mutate(row_id = row_number()) %>%
        filter(row_id <= GOtopN) %>%
        select(-row_id)
    KEGG <- KEGG %>%
        mutate(ONTOLOGY = "KEGG") %>%
        arrange(p.adjust, -Count) %>%
        mutate(row_id = row_number()) %>%
        filter(row_id <= KEGGtopN) %>%
        select(-row_id)
    use_pathway <- rbind(GO, KEGG) %>%
        mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
        dplyr::arrange(ONTOLOGY, desc(p.adjust)) %>%
        mutate(Description = paste0(" ", Description)) %>% 
        mutate(Description = factor(Description, levels = Description)) %>%
        tibble::rowid_to_column("index")

    xaxis_max <- max(-log10(use_pathway$pvalue)) + 1
    pal <- c(KEGG = "#c3e1e6", MF = "#f3dfb7", CC = "#dcc6dc", BP = "#96c38e")

    rect.data <- group_by(use_pathway, ONTOLOGY) %>%
        summarize(n = n()) %>%
        ungroup() %>%
        mutate(
            ymax = cumsum(n),
            ymin = lag(ymax, default = 0) + 0.6,
            ymax = ymax + 0.4
        )
    
    p <- use_pathway %>%
        ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
        geom_round_col(
            aes(y = Description),
            width = 0.6, alpha = 1 # 这里改了一下透明度 原来是0.8
        ) +
        geom_text(
            aes(x = 0.05, label = paste0(" ", Description)),
            # aes(x = 0.05, label = Description),
            hjust = 0, size = 5
        ) +
        geom_point(
            aes(x = -0.02 * xaxis_max, size = Count),
            shape = 21,
            show.legend = F
        ) +
        geom_text(
            aes(x = -0.02 * xaxis_max, label = Count)
        ) +
        scale_size_continuous(name = "Count", range = c(5, 12)) +
        geom_segment(
            aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
            data = data.frame(xaxis_max = xaxis_max),
            linewidth = 1.5,
            inherit.aes = FALSE
        ) +
        # labs(y = NULL) +
        scale_x_continuous(name = expression(-log[10](pvalue))) + # 设置 x 轴标题
        scale_fill_manual(name = "Category", values = pal, drop = F) +
        # scale_colour_manual(values = pal, drop = F) +
        guides(fill = guide_legend(reverse = TRUE)) +
        theme_prism() +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_text(),
            axis.text.x = element_text(size = 14),
            axis.title.x = element_text(size = 14)
        )

    print(p)
    ggsave("/share/home/lik/Scripts/likun/Random_tasks/Results/WGCNA/ace/log_impute_scale_INT/median/do_not_redistribution/keggo_plots_dfs/lmh.intensity_pink_beautiful_keggo.pdf", p, width = 7, height = 18)
}

module_highlight <- function(go_enrich = NULL, kegg_enrich = NULL, GOtopN = 15, KEGGtopN = 15, useOn = NULL, n_space = 0, title = NULL) {
    pal <- c(KEGG = "#c3e1e6", MF = "#f3dfb7", CC = "#dcc6dc", BP = "#96c38e")
    # library(gground)  这个包好像不需要, 尴尬了
    library(ggprism)
    GO <- NULL
    KEGG <- NULL
    if (!is.null(go_enrich)) {
      if (class(go_enrich)[[1]] == "enrichResult") {
        ego_readable <- setReadable(go_enrich, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
        GO <- as.data.frame(ego_readable)
      } else {
        GO <- go_enrich
      }
        if (nrow(GO) > 0) {
            GO$Description <- paste0(" ", GO$Description)
        }
    }
    if (!is.null(kegg_enrich)) {
      if (class(kegg_enrich)[[1]] == "enrichResult") {
        ekegg_readable <- setReadable(kegg_enrich, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
        KEGG <- as.data.frame(ekegg_readable)
      } else {
        KEGG <- kegg_enrich
      }
        if (nrow(KEGG) > 0) {
            KEGG$Description <- paste0(" ", KEGG$Description)
        }
    }
    if ((!is.null(GO)) & (!is.null(KEGG))) {
        GO <- GO %>%
            group_by(ONTOLOGY) %>%
            arrange(p.adjust, -Count) %>%
            mutate(row_id = row_number()) %>%
            filter(row_id <= GOtopN) %>%
            select(-row_id)
        KEGG <- KEGG %>%
            mutate(ONTOLOGY = "KEGG") %>%
            arrange(p.adjust, -Count) %>%
            mutate(row_id = row_number()) %>%
            filter(row_id <= KEGGtopN) %>%
            select(-row_id)
        use_pathway <- rbind(GO, KEGG) %>%
            mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
            dplyr::arrange(ONTOLOGY, desc(p.adjust)) %>%
            mutate(Description = factor(Description, levels = Description)) %>%
            tibble::rowid_to_column("index")
        # use_pathway <- GO %>%
        #     group_by(p.adjust) %>%
        #     top_n(1, wt = Count) %>%
        #     group_by(ONTOLOGY) %>%
        #     top_n(GOtopN, wt = -p.adjust) %>%
        #     rbind(
        #         top_n(KEGG, KEGGtopN, -p.adjust) %>%
        #             mutate(ONTOLOGY = "KEGG")
        #     ) %>%
        #     ungroup() %>%
        #     mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
        #     dplyr::arrange(ONTOLOGY, desc(p.adjust)) %>%
        #     mutate(Description = factor(Description, levels = Description)) %>%
        #     tibble::rowid_to_column("index")
    } else if (!is.null(GO)) {
        use_pathway <- GO %>%
            group_by(p.adjust) %>%
            top_n(1, wt = Count) %>%
            group_by(ONTOLOGY) %>%
            top_n(GOtopN, wt = -p.adjust) %>%
            ungroup() %>%
            mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
            dplyr::arrange(ONTOLOGY, desc(p.adjust)) %>%
            mutate(Description = factor(Description, levels = Description)) %>%
            tibble::rowid_to_column("index")
    } else if (!is.null(KEGG)) {
        use_pathway <- KEGG %>%
            mutate(ONTOLOGY = "KEGG") %>%
            group_by(ONTOLOGY) %>%
            top_n(KEGGtopN, wt = -p.adjust) %>%
            ungroup() %>%
            mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF", "KEGG")))) %>%
            dplyr::arrange(ONTOLOGY, desc(p.adjust)) %>%
            mutate(Description = factor(Description, levels = Description)) %>%
            tibble::rowid_to_column("index")
    } else {
        message("Both GO and KEGG enrichment results are NULL.")
        return(NULL)
    }
    if (!is.null(useOn)) {
        use_pathway <- use_pathway %>%
            filter(ONTOLOGY %in% useOn) %>%
            dplyr::select(-index) %>%
            rowid_to_column("index")
    }
    if (nrow(use_pathway) == 0) {
        message("No significant pathways found.")
        return(NULL)
    }
    xaxis_max <- max(-log10(use_pathway$pvalue)) + 1

    rect.data <- group_by(use_pathway, ONTOLOGY) %>%
        summarize(n = n()) %>%
        ungroup() %>%
        mutate(
            ymax = cumsum(n),
            ymin = lag(ymax, default = 0) + 0.6,
            ymax = ymax + 0.4
        )
    p <- use_pathway %>%
        ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
        geom_col(
            aes(y = Description),
            width = 0.6, alpha = 1 # 这里改了一下透明度 原来是0.8
        ) +
        geom_text(
            aes(x = 0.05, label = paste0(" ", Description)),
            hjust = 0, size = 5
        ) +
        geom_point(
            aes(x = -0.02 * xaxis_max, size = Count),
            shape = 21,
            show.legend = F
        ) +
        geom_text(
            aes(x = -0.02 * xaxis_max, label = Count)
        ) +
        scale_size_continuous(name = "Count", range = c(5, 12)) +
        geom_segment(
            aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
            data = data.frame(xaxis_max = xaxis_max),
            linewidth = 1.5,
            inherit.aes = FALSE
        ) +
        # labs(y = NULL) +
        scale_x_continuous(name = expression(-log[10](pvalue))) + # 设置 x 轴标题
        scale_fill_manual(name = "Category", values = pal, drop = F) +
        scale_colour_manual(values = pal, drop = F) +
        guides(fill = guide_legend(reverse = TRUE)) +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.x = element_line(),
            axis.line = element_blank(),
            axis.ticks.y = element_blank(),
            legend.title = element_text(),
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.title.x = element_text(size = 14)
        )
    if (!is.null(title)){
      p <- p + labs(title = title) + theme(plot.title = element_text(hjust = 0.5, size = 14))
    }
    print(p)

    return(p = p)
}






my_boxplot <- function(df, var_name = c("ADNC_LMH", "Braak_NFT_stage"), show = c("medians", "totals"), scntfc_scle = T, range = NULL, title = NULL) {
    # 输入一个数据框，列名是medians或者totals以及47个协变量，是用来画medians或者totals随着ADNC_LMH或者Braak_NFT_stage变化的箱线图
    # 参数匹配
    var_name <- match.arg(var_name)
    show <- match.arg(show)
    # 复制数据框以避免修改原始数据
    x <- df
    # 将 ADNC_LMH 或其他变量转换为因子
    if (var_name == "ADNC_LMH") {
        x[[var_name]] <- factor(x[[var_name]],
            levels = c(0, 1, 2, 3),
            labels = c("HC", "ADNC_L", "ADNC_M", "ADNC_H")
        )
    } else if (var_name == "Braak_NFT_stage") {
        # 假设 Braak_NFT_stage 也有类似的分类，需根据实际数据调整
        x[[var_name]] <- factor(x[[var_name]], levels = c(0, 1, 2, 3, 4, 5, 6))
    }
    # 绘制箱线图
    if (scntfc_scle == T) {
        p <- x %>%
            ggplot(aes(x = .data[[var_name]], y = .data[[show]])) +
            geom_boxplot(aes(fill = .data[[var_name]]), width = 0.6, position = position_dodge(width = 0.1), outlier.shape = NA) +
            geom_jitter(color = "black", alpha = 0.3, width = 0.2, size = 1.1) +
            scale_y_continuous(labels = scales::scientific) + # 添加科学计数法
            theme_bw() +
            labs(x = var_name, y = show) +
            ggtitle(label = title) +
            theme(
                panel.grid.major = element_line(linewidth = 0.5, linetype = "solid", color = "lightgray"),
                panel.grid.minor = element_line(linewidth = 0.25, linetype = "dashed", color = "lightgray"),
                panel.grid.major.x = element_line(linewidth = 0.5),
                panel.grid.major.y = element_line(linewidth = 0.5),
                panel.grid.minor.x = element_line(linewidth = 0.5),
                panel.grid.minor.y = element_line(linewidth = 0.5),
                legend.position = "none",
                axis.text = element_text(size = 25),
                axis.title = element_text(size = 25),
                plot.title = element_text(size = 25, hjust = 0.5)
            )
    } else {
        p <- x %>%
            ggplot(aes(x = .data[[var_name]], y = .data[[show]])) +
            geom_boxplot(aes(fill = .data[[var_name]]), width = 0.6, position = position_dodge(width = 0.1), outlier.shape = NA) +
            geom_jitter(color = "black", alpha = 0.3, width = 0.2, size = 1.1) +
            # scale_y_continuous(labels = scales::scientific) +  # 添加科学计数法
            theme_bw() +
            labs(x = var_name, y = show) +
            ggtitle(label = title) +
            theme(
                panel.grid.major = element_line(linewidth = 0.5, linetype = "solid", color = "lightgray"),
                panel.grid.minor = element_line(linewidth = 0.25, linetype = "dashed", color = "lightgray"),
                panel.grid.major.x = element_line(linewidth = 0.5),
                panel.grid.major.y = element_line(linewidth = 0.5),
                panel.grid.minor.x = element_line(linewidth = 0.5),
                panel.grid.minor.y = element_line(linewidth = 0.5),
                legend.position = "none",
                axis.text = element_text(size = 25),
                axis.title = element_text(size = 25),
                plot.title = element_text(size = 25, hjust = 0.5)
            )
    }
    if (!is.null(range)) {
        p <- p + coord_cartesian(ylim = range)
    }
    # 返回图形对象
    return(p)
}





# library(ggplot2)
# library(cowplot)


# library(trackViewer)
# features <- GRanges("chr1", IRanges(c(1, 501, 1001),
#                                     width=c(120, 400, 405),
#                                     names=paste0("block", 1:3)),
#                     fill = c("#FF8833", "#51C6E6", "#DFA32D"),
#                     height = c(0.02, 0.05, 0.08))
# SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
# sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)),
#                      color = sample.int(6, length(SNP), replace=TRUE),
#                      score = rep(0, length(SNP)))
# lolliplot(sample.gr, features)
