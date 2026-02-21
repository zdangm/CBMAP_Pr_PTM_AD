suppressPackageStartupMessages({
  library(data.table)
  library(writexl)
  library(openxlsx)
  library(magrittr)
  library(readxl)
  library(patchwork)
  library(limma)
  library(UpSetR)
  library(gridExtra)
  library(qs2)
  library(optparse)
  library(purrr)
  library(ggrepel)
  library(impute)
  library(clusterProfiler)
  library(tidyverse)
  library(conflicted)
  library(ggpubr)
  options(tibble.width = Inf)
})
rm(list = ls())
conflict_prefer_all("dplyr")

if (F) {
  dm <- read_tsv("domains.tab", skip = 3) %>% as.data.frame()
  dm <- rbind(colnames(dm), dm)
  colnames(dm) <- "tmp"
  res <- data.frame()
  for (i in seq.int(nrow(dm))) {
    items <- str_split_1(dm[i, ], " ")
    items <- items[items != ""]
    res <- rbind(res, items)
    message(i)
  }
  res <- res[, 1:22]
  colname <- str_split_1(
    "target_name        accession1   tlen query_name           accession2   qlen   seq_E-value  seq_score  seq_bias   dom_No  dom_of  dom_c-Evalue  dom_i-Evalue  dom_score  dom_bias  hmm_from    hmm_to  ali_from    ali_to  env_from    env_to  acc description_of_target",
    " "
  )
  colname <- colname[colname != ""]
  colnames(res) <- colname[1:22]
  res <- as.tibble(res)
  n_distinct(res$target_name)
  n_distinct(res$query_name)

  # dup_res <- res %>% distinct(target_name, query_name)
  res <- res %>% mutate(accPF = str_split_i(accession2, "\\.", 1))
  res %>% filter(query_name == "Filament_head")
  res %>% filter(target_name == "P10636")

  clan <- read_tsv("../Pfam_hmmer/PfamA/Pfam-A.clans.tsv", col_names = F) %>%
    rename(
      accPF = X1,
      accCl = X2,
      nameCl = X3,
      namePF = X4,
      Description = X5
    )
  res <- res %>% left_join(clan, by = "accPF")
  write_csv(res, file = "../processed/domain2clan.csv")
}


clan <- read_tsv("../Pfam_hmmer/PfamA/Pfam-A.clans.tsv", col_names = F) %>%
  rename(
    accPF = X1,
    accCl = X2,
    nameCl = X3,
    namePF = X4,
    Description = X5
  )
res <- read_csv("../processed/domain2clan.csv")
res <- res %>%
  select(
    -accession1,
    -tlen,
    -accession2,
    -(qlen:hmm_to),
    -acc,
    -query_name
  )
res %<>%
  # select(-query_name) %>%
  rename(
    accPr = target_name
  )
cbmap_ph <- read_excel("result_table.xlsx", sheet = 20, skip = 1) %>%
  select(id, ADNC.LMH.num.adj.P.Val) %>%
  mutate(
    accPr = str_split_i(str_split_i(id, "_", 1), "\\.", 2),
    pos = str_split_i(id, "_", 2)
  )
cbmap_ub <- read_excel("result_table.xlsx", sheet = 26, skip = 1) %>%
  select(id, ADNC.LMH.num.adj.P.Val) %>%
  mutate(
    accPr = str_split_i(str_split_i(id, "_", 1), "\\.", 2),
    pos = str_split_i(id, "_", 2)
  )
cbmap_ac <- read_excel("result_table.xlsx", sheet = 29, skip = 1) %>%
  select(id, ADNC.LMH.num.adj.P.Val) %>%
  mutate(
    accPr = str_split_i(str_split_i(id, "_", 1), "\\.", 2),
    pos = str_split_i(id, "_", 2)
  )
cbmap.ls <- list(
  ph = cbmap_ph,
  ub = cbmap_ub,
  ac = cbmap_ac
)

res <- res %>%
  mutate(
    ali_from = as.numeric(ali_from),
    ali_to = as.numeric(ali_to),
    env_from = as.numeric(env_from),
    env_to = as.numeric(env_to)
  ) %>%
  filter(!is.na(ali_from), !is.na(ali_to), !is.na(env_from), !is.na(env_to))


final.ls <- list()
for (i in seq_along(cbmap.ls)) {
  cbmap_res <- cbmap.ls[[i]]
  cbmap_res$pos <- as.numeric(cbmap_res$pos)
  tmp <- data.frame()
  for (ii in seq.int(nrow(cbmap_res))) {
    this_accpr <- cbmap_res$accPr[ii]
    this_pos <- cbmap_res$pos[ii]

    sub_res <- res %>%
      filter(accPr == this_accpr, this_pos <= env_to, this_pos >= env_from)
    if (nrow(sub_res) == 0) {
      tmp <- bind_rows(
        tmp,
        bind_cols(
          cbmap_res[ii, ],
          ali_from = NA,
          ali_to = NA,
          env_from = NA,
          env_to = NA,
          accPF = NA,
          accCl = NA,
          nameCl = NA,
          namePF = NA,
          Description = NA
        )
      )
    } else if (nrow(sub_res) > 1) {
      print(cbmap_res$id[ii])
      print(sub_res)
      sub_res <- sub_res[1, ] %>% select(-accPr)
      tmp <- bind_rows(tmp, bind_cols(cbmap_res[ii, ], sub_res))
    } else {
      sub_res <- sub_res %>% select(-accPr)
      tmp <- bind_rows(tmp, bind_cols(cbmap_res[ii, ], sub_res))
    }
  }
  final.ls[[i]] <- tmp
}
names(final.ls) <- names(cbmap.ls)

sub.fin.ls <- list()
for (i in seq_along(final.ls)) {
  tmp <- final.ls[[i]] %>%
    filter(!is.na(namePF)) %>%
    # filter(!is.na(nameCl)) %>%
    as_tibble() %>%
    mutate(
      is_DEP = if_else(ADNC.LMH.num.adj.P.Val <= 0.05, 1, 0)
    ) %>%
    select(-ADNC.LMH.num.adj.P.Val)
  sub.fin.ls[[i]] <- tmp
}
names(sub.fin.ls) <- names(final.ls)
sub.fin.ls[[1]]
lapply(sub.fin.ls, dim)
# [1] 3507   13
# [1] 9160   13
# [1] 3337   13
sub.fin.ls[[1]]$nameCl %>% n_distinct()
# [1] 280
sub.fin.ls[[1]]$namePF %>% n_distinct()
# [1] 884
sub.fin.ls[[2]]$nameCl %>% n_distinct()
# [1] 336
sub.fin.ls[[2]]$namePF %>% n_distinct()
# [1] 1386
sub.fin.ls[[3]]$nameCl %>% n_distinct()
# [1] 270
sub.fin.ls[[3]]$namePF %>% n_distinct()
# [1] 825

all_or_not <- "all"
# all_or_not <- "separate"

if (all_or_not == "all") {
  df.fin <- bind_rows(sub.fin.ls)
  df.fin <- df.fin %>%
    mutate(Source = "all")
} else if (all_or_not == "separate") {
  sub.fin.ls.bak <- sub.fin.ls
  for (i in seq_along(sub.fin.ls)) {
    ptm <- names(sub.fin.ls)[i]
    x <- sub.fin.ls[[i]]
    x <- x %>% mutate(Source = ptm)
    sub.fin.ls[[i]] <- x
  }

  df.fin <- bind_rows(sub.fin.ls)
}
df.fin.bak <- df.fin


test.df <- tibble(
  PF = character(),
  Fisher_P = numeric(),
  Fisher_beta = numeric(),
  mat = list(),
  Source = character()
)
PFs <- unique(df.fin$accPF)
Sources <- unique(df.fin$Source)
for (ii in seq_along(PFs)) {
  for (j in seq_along(Sources)) {
    l_u <- df.fin %>%
      filter(accPF == PFs[ii], is_DEP == 1, Source == Sources[j]) %>%
      nrow()
    r_u <- df.fin %>%
      filter(accPF != PFs[ii], is_DEP == 1, Source == Sources[j]) %>%
      nrow()
    l_l <- df.fin %>%
      filter(accPF == PFs[ii], is_DEP == 0, Source == Sources[j]) %>%
      nrow()
    r_l <- df.fin %>%
      filter(accPF != PFs[ii], is_DEP == 0, Source == Sources[j]) %>%
      nrow()
    mat <- matrix(
      data = c(l_u, r_u, l_l, r_l),
      nrow = 2,
      byrow = T,
      dimnames = list(c("DEP", "not_DEP"), c("in_PF", "not_in"))
    )
    ft <- fisher.test(mat)
    test.df <- test.df %>%
      bind_rows(
        tibble(
          PF = PFs[ii],
          Fisher_P = ft$p.value,
          Fisher_beta = ft$estimate,
          mat = list(mat),
          Source = Sources[j]
        )
      )
  }
}

test.df <- test.df %>%
  left_join(clan, by = c("PF" = "accPF")) %>%
  # left_join(clan, by = c("PF" = "accCl")) %>%
  arrange(Fisher_P)

test.df <- test.df %>%
  group_by(Source) %>%
  mutate(Fisher_fdr = p.adjust(Fisher_P)) %>%
  ungroup()
test.df %>% filter(Fisher_fdr < 0.05)

sig_PFs <- test.df %>% filter(Fisher_fdr < 0.05) %>% pull(PF) %>% unique()
df.fin %>% filter(is_DEP == 1, accPF == sig_PFs[3])


fwrite(test.df, file = "test.df.txt")

mat.ls <- test.df[test.df$Fisher_fdr < 0.05, ]
for (i in seq.int(nrow(mat.ls))) {
  des <- mat.ls$Description[i]
  mat <- mat.ls[i, ][["mat"]][[1]]
  ptm <- mat.ls[i, ]$Source
  message(paste0("Description: ", des, ". In: ", ptm))
  print(mat)
}


test.df <- fread("test.df.txt")

test.df %>% filter(Fisher_fdr < 0.05) %>% pull(PF)

sig_PF_df <-
  data.frame(
    PF = test.df %>% filter(Fisher_fdr < 0.05) %>% arrange(Source) %>% pull(PF),
    source = c(rep("ac", 2), rep("ph", 1), rep("ub", 12))
  )

res <- data.frame()
for (i in seq.int(nrow(sig_PF_df))) {
  pf <- sig_PF_df[i, 1]
  src <- sig_PF_df[i, 2]
  tmp <- df.fin %>%
    filter(accPF == pf, is_DEP == 1, Source == src)
  res <- rbind(res, tmp)
}
res
res %>% write_xlsx("sig_sites_in_sig_PFs_groupby_ptm.xlsx")

sig.test.df <- test.df %>%
  filter(Fisher_fdr < 0.05) %>%
  rename(
    ONTOLOGY = Source,
    p.adjust = Fisher_fdr,
    pvalue = Fisher_P
  )
go_enrich <- sig.test.df

go_enrich$Description %>% table()
which(go_enrich$Description == "Tau and MAP protein, tubulin-binding repeat")
go_enrich$Description[5] <- "Tau and MAP protein, tubulin-binding repeat "

pal <- c(ph = "#fbb4ae", ub = "#b3cde3", ac = "#ccebc5")
library(ggprism)
GO <- NULL
if (!is.null(go_enrich)) {
  GO <- go_enrich
  # if (nrow(GO) > 0) {
  #   GO$Description <- paste0(" ", GO$Description)
  # }
}


rect.data <- group_by(use_pathway, ONTOLOGY) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )


use_pathway <- GO %>%
  mutate(
    ONTOLOGY = factor(ONTOLOGY, levels = c("ac", "ub", "ph"), ordered = T)
  ) %>%
  dplyr::arrange(ONTOLOGY, desc(p.adjust)) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column("index")
xaxis_max <- max(-log10(use_pathway$pvalue)) + 1


p <- use_pathway %>%
  ggplot(aes(-log10(pvalue), y = index, fill = ONTOLOGY)) +
  geom_col(
    aes(y = Description),
    width = 0.6,
    alpha = 1 #
  ) +
  geom_text(
    aes(x = 0.05, label = paste0(" ", Description)),
    hjust = 0,
    size = 5
  ) +
  geom_segment(
    aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
    data = data.frame(xaxis_max = xaxis_max),
    linewidth = 1.5,
    inherit.aes = FALSE
  ) +
  scale_x_continuous(name = expression(-log[10](P))) + #
  scale_fill_manual(name = "Category", values = pal, drop = F) +
  scale_colour_manual(values = pal, drop = F) +
  guides(fill = guide_legend(reverse = T)) +
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

p <- p +
  labs(title = "Enrichment in Protein Family") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

print(p)
ggsave(filename = "Enrichment_in_PF.pdf", p, width = 8, height = 6)
