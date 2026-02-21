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
  library(ggrepel)
  library(impute)
  library(readxl)
  library(clusterProfiler)
  library(tidyverse)
  library(conflicted)
  library(mediation)
  library(sva)
  options(tibble.width = Inf)
})
rm(list = ls())
conflict_prefer_all("dplyr")


RES_SUB_DIR <- "../Results/ptm_cor_prWGCNA/"

sly <- read_rds("adnc_limma_all_586.rds")
MEs <- sly$MEs
dim(MEs)
MEs <- MEs %>%
  select(ME4, ME11, ME16)
id_ord <- rownames(MEs)


source("helper_func.R")
sig.sites.ls <- list(
  phospho = "phos/All.sites.All.traits.intensity.limma.xlsx",
  ubiq = "ubiq/All.sites.All.traits.intensity.limma.xlsx",
  ace = "ace/All.sites.All.traits.intensity.limma.xlsx"
)
intensity.ls <- list(
  phospho = "phos/mis_log2_impt_sc_rint.729.all.info.ls.QS2",
  ubiq = "ubiq/mis_log2_impt_sc_rint.729.all.info.ls.QS2",
  ace = "ace/mis_log2_impt_sc_rint.729.all.info.ls.111.QS2"
)


option_list <- list(
  make_option(
    c("--n_batch"),
    type = "integer",
    default = 10,
    help = "total",
    metavar = "number"
  ),
  make_option(
    c("--batch"),
    type = "integer",
    default = 1,
    help = "index",
    metavar = "number"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


batch_id <- opt$batch
n_batch <- opt$n_batch


MEs <- MEs %>%
  rownames_to_column("id") %>%
  filter(id != "PTB466") %>%
  column_to_rownames("id")
M.ls <- list(
  phospho = NULL,
  ubiq = NULL,
  ace = NULL
)
# ace -----------------------------------------------------------------
sig.sites <- read_excel(sig.sites.ls[[3]]) %>%
  filter(ADNC.LMH.num.adj.P.Val < 0.05) %>%
  pull(id)
qs_readm(intensity.ls[[3]])

M.ls$ace <- log1p_raw.all.info.ls[[1]] %>%
  select(all_of(sig.sites), id, sex_male, RIN, age, PMD, ADNC_LMH) %>%
  filter(id %in% rownames(MEs)) %>%
  arrange(match(id, rownames(MEs))) %>%
  bind_cols(MEs)

a.result_tbl <- expand_grid(X = colnames(MEs), PTM = "ace", M = sig.sites) %>%
  mutate(result = vector("list", n()))


# phospho -----------------------------------------------------------------
sig.sites <- read_excel(sig.sites.ls[[1]]) %>%
  filter(ADNC.LMH.num.adj.P.Val < 0.05) %>%
  pull(id)
qs_readm(intensity.ls[[1]])

M.ls$phospho <- log1p_raw.all.info.ls[[1]] %>%
  select(all_of(sig.sites), id, sex_male, RIN, age, PMD, ADNC_LMH) %>%
  filter(id %in% rownames(MEs)) %>%
  arrange(match(id, rownames(MEs))) %>%
  bind_cols(MEs)

p.result_tbl <- expand_grid(
  X = colnames(MEs),
  PTM = "phospho",
  M = sig.sites
) %>%
  mutate(result = vector("list", n()))


# ubiq -----------------------------------------------------------------
sig.sites <- read_excel(sig.sites.ls[[2]]) %>%
  filter(ADNC.LMH.num.adj.P.Val < 0.05) %>%
  pull(id)
qs_readm(intensity.ls[[2]])

M.ls$ubiq <- log1p_raw.all.info.ls[[1]] %>%
  select(all_of(sig.sites), id, sex_male, RIN, age, PMD, ADNC_LMH) %>%
  filter(id %in% rownames(MEs)) %>%
  arrange(match(id, rownames(MEs))) %>%
  bind_cols(MEs)

u.result_tbl <- expand_grid(X = colnames(MEs), PTM = "ubiq", M = sig.sites) %>%
  mutate(result = vector("list", n()))


# bind --------------------------------------------------------------------

result_tbl.bak <- rbind(p.result_tbl, u.result_tbl, a.result_tbl)
result_tbl <- result_tbl.bak %>% mutate(group = ntile(row_number(), n_batch))
result_tbl <- result_tbl %>% filter(group == batch_id)


# Mediation now -----------------------------------------------------------

for (k in seq.int(nrow(result_tbl))) {
  ptm <- result_tbl$PTM[k]
  ptm_site <- result_tbl$M[k]
  me <- result_tbl$X[k]

  meddata <- M.ls[[ptm]] %>%
    select(ADNC_LMH, all_of(me), all_of(ptm_site), sex_male:PMD)
  colnames(meddata)[1:3] <- c("Y", "X", "M")
  if (ncol(meddata) != 7) {
    stop(sprintf(
      "colnames of meddata are: %s",
      paste(colnames(meddata), collapse = ", ")
    ))
  }

  fitYX <- lm(Y ~ X + sex_male + age + PMD + RIN, data = meddata)
  fitYM <- lm(Y ~ M + sex_male + age + PMD + RIN, data = meddata)

  fitY <- lm(Y ~ X + M + sex_male + age + PMD + RIN, data = meddata)
  fitM <- lm(M ~ X + sex_male + age + PMD + RIN, data = meddata)

  res_raw <- as.data.frame(summary(fitYX)$coefficients)[2, ]
  res_reg <- as.data.frame(summary(fitY)$coefficients)[2, ]

  res_YM <- as.data.frame(summary(fitYM)$coefficients)[2, ]
  res_MX <- as.data.frame(summary(fitM)$coefficients)[2, ]

  # tic("start")
  med_result <- mediate(
    model.m = fitM,
    model.y = fitY,
    treat = "X",
    mediator = "M",
    covariates = meddata[, c("sex_male", "age", "PMD", "RIN")],
    boot = TRUE,
    sims = 1000
    # sims = 50000
  )
  # toc()

  res <- data.frame(
    Effect = c("ACME", "ADE", "Total Effect", "Prop. Mediated", "M_X", "Y_M"),
    Estimate = c(
      med_result$d0, # ACME
      med_result$z0, # ADE   res_reg$Estimate,
      med_result$tau.coef, # Total Effect  res_raw$Estimate,
      med_result$n0, # Prop. Mediated
      res_MX$Estimate,
      res_YM$Estimate
    ),
    SE = c(
      NA,
      res_reg$`Std. Error`,
      res_raw$`Std. Error`,
      NA,
      res_MX$`Std. Error`,
      res_YM$`Std. Error`
    ),
    CI_Lower = c(
      med_result$d0.ci[1], # ACME 95% CI
      med_result$z0.ci[1], # ADE 95% CI
      med_result$tau.ci[1], # Total 95% CI
      med_result$n0.ci[1], # Prop. Mediated 95% CI
      NA,
      NA
    ),
    CI_Upper = c(
      med_result$d0.ci[2], # ACME 95% CI
      med_result$z0.ci[2], # ADE 95% CI
      med_result$tau.ci[2], # Total 95% CI
      med_result$n0.ci[2], # Prop. Mediated 95% CI
      NA,
      NA
    ),
    p_value = c(
      med_result$d0.p, # ACME p-value
      res_reg$`Pr(>|t|)`, # ADE p-value  med_result$z0.p
      res_raw$`Pr(>|t|)`, # Total p-value
      med_result$n0.p, # Prop. Mediated p-value
      res_MX$`Pr(>|t|)`,
      res_YM$`Pr(>|t|)`
    )
  )

  row_idx <- result_tbl %>%
    mutate(row_id = row_number()) %>%
    filter(X == me, M == ptm_site, PTM == ptm) %>%
    pull(row_id)
  result_tbl$result[[row_idx]] <- res
}


qs_save(
  result_tbl,
  file = paste0(RES_SUB_DIR, "batch", batch_id, "_of_", n_batch, ".qs")
)


if (F) {
  batches <- list.files(path = RES_SUB_DIR, pattern = "qs", full.names = T)
  final_df <- data.frame()
  for (i in batches) {
    df <- qs_read(i)
    message(i)
    for (j in seq.int(nrow(df))) {
      pairship <- df[j, ] %>% unnest(cols = c(result))
      pairship$Effect <- c("ACME", "ADE", "TTLE", "PropM", "MX", "YM")

      pairship <- pairship %>%
        select(-group) %>%
        pivot_wider(
          id_cols = X:M,
          names_from = Effect,
          values_from = Estimate:p_value
        )

      final_df <- final_df %>% bind_rows(pairship)
    }
  }

  write_xlsx(final_df, path = paste0(RES_SUB_DIR, "ptm_med_ME_LMH.xlsx"))
}


if (F) {
  final_df <- read_excel(path = paste0(RES_SUB_DIR, "ptm_med_ME_LMH.xlsx"))

  final_df %>% filter(p_value_YM < 0.05) %>% nrow()
  dim(final_df)

  {
    final_df <- final_df %>%
      filter(
        fdr_ACME < 0.05,
        # fdr_TTLE < 0.05,
        # fdr_YM < 0.05,
        fdr_MX < 0.05
      )
    dim(final_df)
    tmp <- final_df %>%
      group_by(X, PTM) %>%
      summarise(
        n_sites = n(),
        genes = list(str_split_i(M, pattern = "\\.", i = 1))
      ) %>%
      arrange(PTM, X)

    vec <- sapply(tmp$genes, function(x) {
      x <- unique(x)
      length(x)
    })
    tmp %>% bind_cols(n_genes = vec) %>% select(-genes)
  }

  for_enrich <- final_df %>%
    group_by(
      X,
      PTM
    ) %>%
    summarise(
      nrow = n(),
      genes = list(str_split_i(M, pattern = "\\.", i = 1))
    )
  i <- 1

  wb <- createWorkbook()
  wb1 <- createWorkbook()
  for (i in seq.int(nrow(for_enrich))) {
    pairship <- paste0(for_enrich[i, ]$X, "_", for_enrich[i, ]$PTM)

    {
      tmp <- final_df %>%
        filter(
          X == for_enrich[i, ]$X,
          PTM == for_enrich[i, ]$PTM,
        ) %>%
        select(
          X:M,
          everything(),
          -starts_with("CI"),
          -starts_with("SE"),
          -ends_with("TTLE")
        ) %>%
        arrange(
          p_value_ACME,
          p_value_MX
        )
      tmp$loc_on_tau <- sapply(tmp$M, from_longest_isof_to_clinical_isof)
      tmp <- tmp %>%
        mutate(M = if_else(is.na(loc_on_tau), M, loc_on_tau)) %>%
        select(-loc_on_tau)

      addWorksheet(wb1, sheetName = pairship)
      writeData(wb1, sheet = pairship, tmp)
    }

    {
      sig_genes <- for_enrich[i, ]$genes[[1]] %>% unique()

      go.res <- enrichment_gene(
        sig_genes,
        var_name = "",
        enrich_method = "GO",
        PTM = for_enrich[i, ]$PTM
      )
      kegg.res <- enrichment_gene(
        sig_genes,
        var_name = "",
        enrich_method = "KEGG",
        PTM = for_enrich[i, ]$PTM
      )
      kegg.res <- as.data.frame(kegg.res)
      pic <- module_highlight(
        go_enrich = go.res,
        kegg_enrich = kegg.res,
        GOtopN = 15,
        KEGGtopN = 15
      ) +
        ggtitle(sprintf(
          "%s mediates %s to AD",
          for_enrich[i, ]$PTM,
          for_enrich[i, ]$X
        ))
      ggsave(
        paste0(
          RES_SUB_DIR,
          for_enrich[i, ]$PTM,
          "_",
          for_enrich[i, ]$X,
          ".png"
        ),
        pic,
        dpi = 900,
        width = 8,
        height = 10
      )

      addWorksheet(
        wb,
        sheetName = paste0(for_enrich[i, ]$PTM, "_", for_enrich[i, ]$X, "_GO")
      )
      writeData(
        wb,
        sheet = paste0(for_enrich[i, ]$PTM, "_", for_enrich[i, ]$X, "_GO"),
        as.data.frame(go.res)
      )

      addWorksheet(
        wb,
        sheetName = paste0(for_enrich[i, ]$PTM, "_", for_enrich[i, ]$X, "_KEGG")
      )
      writeData(
        wb,
        sheet = paste0(for_enrich[i, ]$PTM, "_", for_enrich[i, ]$X, "_KEGG"),
        as.data.frame(kegg.res)
      )
    }
  }
  saveWorkbook(wb, file = paste0(RES_SUB_DIR, "mediation_keggo.xlsx"))
  saveWorkbook(
    wb1,
    file = paste0(RES_SUB_DIR, "sig_mediations.xlsx"),
    overwrite = T
  )
}


sheets <- excel_sheets(paste0(RES_SUB_DIR, "mediation_keggo.xlsx"))

# 查看 sheet 数量
length(sheets)
GO.df <- data.frame()
kegg.df <- data.frame()
for (i in seq_along(sheets)) {
  sn <- sheets[i]
  enrich_type <- str_split_i(sn, "_", 3)
  ME <- str_split_i(sn, "_", 2)
  ptm <- str_split_i(sn, "_", 1)

  df <- read_excel(paste0(RES_SUB_DIR, "mediation_keggo.xlsx"), sheet = sn)
  if (nrow(df) == 0) {
    next
  }
  df <- df %>% mutate(PTM = ptm, X = ME)
  if (enrich_type == "KEGG") {
    df <- kegg_results_fixing(df)
    kegg.df <- rbind(kegg.df, df)
  } else {
    GO.df <- rbind(GO.df, df)
  }
}

write_xlsx(kegg.df, path = paste0(RES_SUB_DIR, "tmp_kegg.xlsx"))
write_xlsx(GO.df, path = paste0(RES_SUB_DIR, "tmp_go.xlsx"))

ub.wgcna.keggo <- list.files(pattern = "^lmh.*xlsx$")

go.df <- data.frame()
kegg.df <- data.frame()

colorlabel <- labels2colors(seq.int(length(ub.wgcna.keggo)))
for (i in seq_along(ub.wgcna.keggo)) {
  keggo <- ub.wgcna.keggo[i]
  color <- str_split_i(keggo, "\\.", 2) %>% str_split_i("_", 2)
  idx <- which(color == colorlabel)

  sheets <- excel_sheets(keggo)
  for (j in seq_along(sheets)) {
    sn <- sheets[j]
    df <- read_xlsx(keggo, sheet = sn)
    if (str_detect(sn, "GO")) {
      df <- df %>%
        mutate(ONTOLOGY = sn, Mod = paste0("Mod", idx), ModColor = color)
      go.df <- rbind(go.df, df)
    } else {
      df <- df %>% mutate(Mod = paste0("Mod", idx), ModColor = color)
      kegg.df <- rbind(kegg.df, df)
    }
  }
}
write_xlsx(go.df, path = "go.xlsx")
write_xlsx(kegg.df, path = "kegg.xlsx")
