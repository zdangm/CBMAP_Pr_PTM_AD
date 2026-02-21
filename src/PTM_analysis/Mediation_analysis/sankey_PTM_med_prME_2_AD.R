suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
  library(magrittr)
  library(patchwork)
  library(gridExtra)
  library(qs2)
  library(readxl)
  library(tidyverse)
  library(conflicted)
  library(ggalluvial)
  library(networkD3)
  options(tibble.width = Inf)
})
rm(list = ls())
conflict_prefer_all("dplyr")
RES_SUB_DIR <- "../Results/ptm_cor_prWGCNA/"


myelin <- c("myelin sheath", "axolemma")


cytoskeleton <- c(
  "actin binding",
  "tubulin binding",
  "axon cytoplasm",
  "cortical actin cytoskeleton",
  "cortical cytoskeleton",
  "intermediate filament organization",
  "intermediate filament cytoskeleton organization",
  "structural constituent of cytoskeleton",
  "structural constituent of postsynapse"
) %>%
  unique()

energy <- c(
  str_split_1(
    "mitochondrial electron transport, ubiquinol to cytochrome c
oxidative phosphorylation
aerobic respiration
cellular respiration
aerobic electron transport chain",
    "\n"
  ),
  str_split_1(
    "mitochondrial respiratory chain complex III
respiratory chain complex III
cytochrome complex
inner mitochondrial membrane protein complex
respiratory chain complex",
    "\n"
  ),
  "2-oxoglutarate metabolic process",
  str_split_1(
    "amino acid catabolic process
organic acid catabolic process
carboxylic acid catabolic process
amino acid metabolic process",
    "\n"
  )
) %>%
  unique()


pr_folding <- c(
  "chaperone-mediated protein folding",
  "protein folding",
  "'de novo' post-translational protein folding",
  "'de novo' protein folding",
  "protein refolding",
  "ATP-dependent protein folding chaperone",
  "unfolded protein binding",
  "protein folding chaperone",
  "heat shock protein binding"
) %>%
  unique()

cate2path <- list(
  myelin = myelin,
  cytoskeleton = cytoskeleton,
  respiration = energy,
  `protein folding` = pr_folding
)


path2gene <- tibble(
  cato = character(),
  # path = character(),
  gene = list()
)


sheetname <- expand_grid(
  ptm = c("phospho", "ubiq", "ace"),
  mod = c("ME4", "ME11", "ME16"),
  path = c("GO")
) %>%
  mutate(
    sheetname = str_c(ptm, mod, path, sep = "_"),
    sheetname2 = str_c(mod, ptm, sep = "_")
  )

for (i in seq_len(nrow(sheetname))) {
  sn <- sheetname$sheetname[i]
  enrich <- read_excel(paste0(RES_SUB_DIR, "mediation_keggo.xlsx"), sheet = sn)
  for (j in seq_along(cate2path)) {
    cate <- names(cate2path)[j]
    path <- cate2path[[j]]
    genes <- enrich %>%
      filter(Description %in% path) %>%
      pull(geneID) %>%
      str_split("/") %>%
      unlist()
    old_row_genes <- path2gene %>%
      filter(cato == cate) %>%
      pull(gene) %>%
      unlist()
    path2gene <- path2gene %>% filter(cato != cate)
    new_row_genes <- unique(c(genes, old_row_genes))
    path2gene <- path2gene %>%
      bind_rows(
        tibble(
          cato = cate,
          # path = path,
          gene = list(new_row_genes)
        )
      )
  }
}

{
  path2gene <- path2gene %>% unnest()
  dim(path2gene)
  rep_genes <- path2gene[path2gene$gene %>% duplicated() %>% which(), ]$gene
  path2gene <- path2gene %>% filter(!(gene %in% rep_genes))
  path2gene <- path2gene %>%
    bind_rows(
      tibble(
        cato = c(
          "cytoskeleton",
          "myelin",
          "myelin",
          "protein folding",
          "cytoskeleton",
          "protein folding"
        ),
        gene = c("PLEC", "ERMN", "MYO1D", "HSP90AA1", "MAPT", "SNCA")
      )
    )
  dim(path2gene)
}


sig_med <- data.frame()
for (i in seq_len(nrow(sheetname))) {
  sn <- sheetname$sheetname2[i]
  sub.sig_med <- read_excel(
    paste0(RES_SUB_DIR, "sig_mediations.xlsx"),
    sheet = sn
  ) %>%
    mutate(gene = str_split_i(M, "\\.", 1)) %>%
    select(X, PTM, M, gene, Estimate_ACME)
  sig_med <- rbind(sig_med, sub.sig_med)
}

final_df <- sig_med %>%
  left_join(path2gene, by = "gene") %>%
  mutate(outcome = "AD", cato = if_else(is.na(cato), "other", cato)) %>%
  select(-M) %>%
  group_by(gene) %>%
  mutate(freq = sum(abs(Estimate_ACME))) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  filter(freq >= 0) %>%
  select(-gene, -Estimate_ACME) %>%
  group_by(
    X,
    PTM,
    cato
  ) %>%
  summarise(
    outcome = "AD",
    freq = round(sum(freq), 2) * 100
  ) %>%
  mutate(
    X = factor(X, levels = c("ME4", "ME11", "ME16"), ordered = T),
    PTM = factor(
      PTM,
      levels = c("phospho", "ubiq", "ace"),
      labels = c("phos", "ubiq", "ace"),
      ordered = T
    ),
    cato = factor(
      cato,
      levels = c(
        "cytoskeleton",
        "myelin",
        "protein folding",
        "respiration",
        "other"
      )
    )
  )

final_df

sankey.pic <- ggplot(
  final_df,
  aes(axis1 = X, axis2 = PTM, axis3 = cato, axis4 = outcome, y = freq)
) +
  geom_alluvium(aes(fill = PTM), width = 1 / 12) +
  scale_fill_manual(
    values = c(
      "ace" = "#c7e2c0",
      "phos" = "#f3b0ab",
      "ubiq" = "#b0c8dd"
    )
  ) +
  geom_stratum(width = 1 / 6, fill = "grey90", color = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(
    limits = c("Module", "PTM", "Category", "Outcome"),
    expand = c(.05, .05)
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  # coord_flip() +
  labs(title = "", y = "Frequency")


sankey.pic
ggsave(
  filename = paste0(RES_SUB_DIR, "sankey_plot.pdf"),
  sankey.pic,
  width = 12,
  height = 8
)


final_df %>%
  group_by(PTM, cato) %>%
  summarise(sum_freq = sum(freq)) %>%
  ungroup() %>%
  group_by(cato) %>%
  mutate(
    ss = sum(sum_freq),
    prpt = sum_freq / sum(sum_freq)
  ) %>%
  arrange(cato)
