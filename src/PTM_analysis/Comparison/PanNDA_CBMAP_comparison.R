suppressPackageStartupMessages({
  library(GSEABase)
  library(data.table)
  library(msigdbr)
  library(openxlsx)
  library(magrittr)
  library(patchwork)
  library(limma)
  library(UpSetR)
  library(gridExtra)
  library(qs2)
  library(purrr)
  library(biomaRt)
  library(ggrepel)
  library(clusterProfiler)
  library(tidyverse)
  library(conflicted)
  library(readxl)
  library(org.Hs.eg.db)
  library(RColorBrewer)
  library(VennDiagram)
  library(ComplexHeatmap)
  library(ComplexUpset)
  options(tibble.width = Inf)
})
rm(list = ls())
gc()
conflict_prefer_all("dplyr")
setwd("/data/projects/lik/CBMAP/Random_tasks/Scripts/")


# ==============================================================================
# ==============================================================================

if (F) {
  pannda.sig.ph <- read_csv("../Other_source_data/AD_phospho.csv", skip = 1) %>%
    separate_longer_delim(ModSites, ",") %>%
    select(-1, -4) %>%
    filter(FDR < 0.05) %>%
    mutate(newid = paste0(Gene, "_", ModSites))
  pannda.sig.ub <- read_csv("../Other_source_data/AD_ub.csv", skip = 1) %>%
    separate_longer_delim(ModSites, ";") %>%
    mutate(
      ModSites = if_else(
        str_detect(ModSites, "^K"),
        ModSites,
        paste0("K", ModSites)
      )
    ) %>%
    select(-1, -4) %>%
    filter(FDR < 0.05) %>%
    mutate(newid = paste0(Gene, "_", ModSites))
}


cbmap.pr <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 2,
  skip = 1
)
# whole proteome

pannda.pr1 <- read_xlsx("../Other_source_data/AD_stats.xlsx", skip = 3) %>%
  mutate(
    ft = str_sub(`Unique identifier`, 1, 2),
    acc = str_split_i(`Unique identifier`, "\\.", 2)
  ) %>%
  filter(
    ft == "sp"
  ) %>%
  rename(id = `Unique identifier`) %>%
  mutate(
    gene_protein = paste0(Gene, "_", acc),
    iso = if_else(
      str_detect(str_split_i(id, "\\.", 3), "^\\d+$"),
      str_split_i(id, "\\.", 3),
      "1"
    )
  ) %>%
  filter(
    iso == "1"
  ) %>%
  select(-iso)
colnames(pannda.pr1) <- c(
  "id",
  "gene",
  "comb_p",
  "fdr",
  "avg_log2FC",
  "avg_log2FC_z",
  "log2FC_1",
  "log2FC_z_1",
  "p_1",
  "log2FC_2",
  "log2FC_z_2",
  "p_2",
  "log2FC_3",
  "log2FC_z_3",
  "p_3",
  "used_in_volcano_plot",
  "ft",
  "acc",
  "gene_protein"
)


# insoluble proteome,
pannda.pr2 <- read_csv("../Other_source_data/AD_insoluble.csv") %>%
  mutate(
    ft = str_sub(Identifier, 1, 2),
    acc = str_split_i(Identifier, "\\.", 2)
  ) %>%
  filter(
    ft == "sp"
  ) %>%
  rename(id = Identifier) %>%
  mutate(
    gene_protein = paste0(Gene, "_", acc),
    iso = if_else(
      str_detect(str_split_i(id, "\\.", 3), "^\\d+$"),
      str_split_i(id, "\\.", 3),
      "1"
    )
  ) %>%
  filter(
    iso == "1"
  ) %>%
  select(-iso)
colnames(pannda.pr2) <- c(
  "id",
  "gene",
  "enrich_factor",
  "p",
  "fdr",
  "log2FC",
  "log2FC_z",
  "ft",
  "acc",
  "gene_protein"
)


literature.pr <- read_csv(
  "../Other_source_data/alzheimer_gene_literature_final.csv"
) %>%
  filter(type != "unknown") %>%
  pull(gene_name)


pannda.sig.pr1.name <- pannda.pr1 %>% filter(fdr < 0.05) %>% pull(gene)
pannda.sig.pr2.name <- pannda.pr2 %>% filter(fdr < 0.05) %>% pull(gene)
cbmap.sig.pr.name <- cbmap.pr %>% filter(adj.P.Val < 0.05) %>% pull(genename)


setdiff(
  cbmap.sig.pr.name,
  unique(c(literature.pr, pannda.sig.pr1.name, pannda.sig.pr2.name))
) %>%
  length()

intersect(
  cbmap.sig.pr.name,
  unique(c(literature.pr, pannda.sig.pr1.name, pannda.sig.pr2.name))
) %>%
  length()


pannda.sig.pr <- rbind(pannda.sig.pr1, pannda.sig.pr2)


cbmap.sig.pr <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 2,
  skip = 1
) %>%
  filter(
    adj.P.Val < 0.05
  )


cbmap.sig.ph <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 20,
  skip = 1
) %>%
  select(
    id,
    ADNC.LMH.num.P.Value,
    ADNC.LMH.num.adj.P.Val,
    `Gene name`,
    ADNC.LMH.num.logFC,
    loc_on_MAPT_8
  ) %>%
  # filter(ADNC.LMH.num.adj.P.Val < 0.05)  %>%
  mutate(
    aa = str_split_i(id, "_", 3),
    sites1 = str_split_i(id, "_", 2),
    sites2 = if_else(
      is.na(loc_on_MAPT_8),
      NA,
      str_split_i(loc_on_MAPT_8, "_", 2)
    )
  ) %>%
  mutate(
    newid = paste0(
      `Gene name`,
      "_",
      if_else(is.na(sites2), paste0(aa, sites1), paste0(aa, sites2))
    )
  )


cbmap.sig.ub <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 26,
  skip = 1
) %>%
  select(
    id,
    ADNC.LMH.num.P.Value,
    ADNC.LMH.num.adj.P.Val,
    `Gene name`,
    ADNC.LMH.num.logFC,
    loc_on_MAPT_8
  ) %>%
  # filter(ADNC.LMH.num.adj.P.Val < 0.05)  %>%
  mutate(
    aa = str_split_i(id, "_", 3),
    sites1 = str_split_i(id, "_", 2),
    sites2 = if_else(
      is.na(loc_on_MAPT_8),
      NA,
      str_split_i(loc_on_MAPT_8, "_", 2)
    )
  ) %>%
  mutate(
    newid = paste0(
      `Gene name`,
      "_",
      if_else(is.na(sites2), paste0(aa, sites1), paste0(aa, sites2))
    )
  )


setdiff(cbmap.sig.ub$newid, pannda.sig.ub$newid) %>% n_distinct()
setdiff(cbmap.sig.ph$newid, pannda.sig.ph$newid) %>% n_distinct()
setdiff(cbmap.sig.pr$gene_protein, pannda.sig.pr1$gene_protein) %>% n_distinct()
setdiff(cbmap.sig.pr$gene_protein, pannda.sig.pr2$gene_protein) %>% n_distinct()
setdiff(cbmap.sig.pr$gene_protein, pannda.sig.pr$gene_protein) %>% n_distinct()


library(ggVennDiagram)

genes_A <- cbmap.sig.pr$gene_protein
genes_B <- pannda.sig.pr1$gene_protein
genes_C <- pannda.sig.pr1 %>%
  filter(FDR < 0.05) %>%
  pull(gene_protein)
venn_list <- list(
  "CBMAP sig pr" = genes_A,
  "PanNDA pr" = genes_B,
  "PanNDA sig pr" = genes_C
)
venn_list <- list(
  "CBMAP sig pr" = genes_A,
  "PanNDA pr\n(insoluble)" = genes_B,
  "PanNDA sig pr\n(insoluble)" = genes_C
)


library(VennDiagram)
library(grid)

v3 <- venn.diagram(
  x = venn_list,
  filename = NULL, #

  fill = c("#E64B35", "#4DBBD5", "#00A087"), #
  alpha = 0.5, #
  col = "white", #

  cex = 1.5, #
  cat.cex = 1.5, #
  cat.fontface = "bold",

  cat.pos = c(-20, 20, 150),
  cat.dist = c(0.05, 0.05, 0.05), #

  main = ""
)

grid.newpage()
grid.draw(v3)


# ==============================================================================
# ============================================================================

setdiff(cbmap.sig.ph$newid, pannda.sig.ph$newid) %>% head()


library(Biostrings)

fasta.raw <- readAAStringSet(
  "../Other_source_data/Homo_sapiens_9606_SP_20231220.fasta"
)
names(fasta.raw) <- str_split_i(names(fasta.raw), " ", 1)

library(clusterProfiler)
library(org.Hs.eg.db)
gene_table <- bitr(
  names(fasta.raw),
  fromType = "UNIPROT",
  toType = "SYMBOL",
  OrgDb = org.Hs.eg.db
)
fasta <- fasta.raw[match(gene_table$UNIPROT, names(fasta.raw))]

identical(names(fasta), gene_table$UNIPROT)


mcols(fasta)$gene_symbol <- gene_table$SYMBOL


pannda.all <- rbind(pannda.sig.ph, pannda.sig.ub)


not_match <- data.frame()
not_match_cnt <- 0
for (i in seq.int(nrow(pannda.all))) {
  aa <- pannda.all[i, ]$ModSites %>% str_extract("[A-Za-z]")
  site <- pannda.all[i, ]$ModSites %>% str_extract("\\d+") %>% as.numeric()
  gene <- pannda.all[i, ]$Gene

  if (!(gene %in% mcols(fasta)$gene_symbol)) {
    # message(gene, " not in fasta!!")
    next
  }
  for (ii in which(mcols(fasta)$gene_symbol == gene)) {
    pr.fasta <- fasta[[ii]] %>% as.character()
    uniprot.aa <- str_sub(pr.fasta, site, site)
    if (aa != uniprot.aa) {
      tmp <- data.frame(gene, site, aa, uniprot.aa)
      not_match <- rbind(not_match, tmp)
      not_match_cnt <- not_match_cnt + 1
    }
  }
}
not_match <- not_match %>%
  as_tibble() %>%
  mutate(newid = paste0(gene, "_", aa, site))

not_match %>% filter(gene != "MAPT", aa == "K")
not_match %>% filter(gene != "MAPT", aa != "K")


# ==============================================================================
# miamiplot pvalue comparison
# ============================================================================

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "asia"
)
ptm <- "pr"
ptm <- "phos"
ptm <- "ubiq"

# pannda_or_meta <- "meta"
pannda_or_meta <- "pannda"


if (ptm == "pr") {
  df1 <- cbmap.sig.pr %>%
    select(genename, P.Value, adj.P.Val) %>%
    rename(
      "Gene name" = genename,
      "ADNC.LMH.num.P.Value" = P.Value,
      "ADNC.LMH.num.adj.P.Val" = adj.P.Val
    ) %>%
    mutate(ADNC.LMH.num.logFC = 1)
  df2 <- NULL
  if (pannda_or_meta == "pannda") {
    df2 <- pannda.sig.pr2 %>%
      mutate(
        `Gene name` = str_split_i(gene_protein, "_", 1),
        ADNC.LMH.num.P.Value = p,
        ADNC.LMH.num.logFC = -1,
        ADNC.LMH.num.adj.P.Val = FDR
      ) %>%
      select(
        `Gene name`,
        ADNC.LMH.num.P.Value,
        ADNC.LMH.num.adj.P.Val,
        ADNC.LMH.num.logFC
      )
  } else if (pannda_or_meta == "meta") {
    df2 <- meta.ad.pr %>%
      select(
        `Gene name` = genename,
        ADNC.LMH.num.P.Value = pval,

        ADNC.LMH.num.adj.P.Val = fdr
      ) %>%
      mutate(
        ADNC.LMH.num.logFC = -1
      )
    df2[1, 1] <- "APP"
  }

  df <- rbind(df1, df2)
  df$identifier <- df$`Gene name`
  df$loc_on_MAPT_8 <- df$`Gene name`
  df$for_volcano <- df$`Gene name`
  genes <- df$`Gene name` %>% unique()
} else if (ptm == "phos") {
  df1 <- cbmap.sig.ph %>%
    mutate(
      ADNC.LMH.num.logFC = 1,
      for_volcano = paste(`Gene name`, sites1, aa, sep = "_")
    ) %>%
    select(
      id,
      loc_on_MAPT_8,
      for_volcano,
      `Gene name`,
      ADNC.LMH.num.logFC,
      ADNC.LMH.num.P.Value,
      ADNC.LMH.num.adj.P.Val
    ) %>%
    rename(
      "identifier" = id,
    )

  df1$loc_on_MAPT_8 <- if_else(
    is.na(df1$loc_on_MAPT_8),
    "NA",
    df1$loc_on_MAPT_8
  )
  df1$loc_on_MAPT_8 <- if_else(
    str_detect(df1$loc_on_MAPT_8, "^not"),
    "NA",
    df1$loc_on_MAPT_8
  )
  df1$for_volcano <- if_else(
    df1$loc_on_MAPT_8 == "NA",
    df1$for_volcano,
    df1$loc_on_MAPT_8
  )
  df1 <- df1 %>% select(-loc_on_MAPT_8)
  df1$for_volcano <- str_replace(df1$for_volcano, "Tau441", "TAU")

  df2 <- pannda.sig.ph %>%
    mutate(
      `Gene name` = Gene,
      ADNC.LMH.num.P.Value = `P-Value`,
      ADNC.LMH.num.logFC = -1,
      ADNC.LMH.num.adj.P.Val = FDR,
      identifier = newid,
      site = stringr::str_extract(ModSites, "\\d+"),
      aa = stringr::str_extract(ModSites, "[A-Za-z]+")
    ) %>%
    mutate(
      for_volcano = paste(`Gene name`, site, aa, sep = "_")
    ) %>%
    select(
      identifier,
      for_volcano,
      `Gene name`,
      ADNC.LMH.num.logFC,
      ADNC.LMH.num.P.Value,
      ADNC.LMH.num.adj.P.Val
    ) %>%
    group_by(identifier) %>%
    filter(ADNC.LMH.num.P.Value == min(ADNC.LMH.num.P.Value)) %>%
    ungroup()
  df2$for_volcano <- str_replace(df2$for_volcano, "MAPT", "TAU")
  df <- rbind(df1, df2)
  genes <- df$`Gene name` %>% unique()
} else if (ptm == "ubiq") {
  df1 <- cbmap.sig.ub %>%
    mutate(
      ADNC.LMH.num.logFC = 1,
      for_volcano = paste(`Gene name`, sites1, aa, sep = "_")
    ) %>%
    select(
      id,
      loc_on_MAPT_8,
      for_volcano,
      `Gene name`,
      ADNC.LMH.num.logFC,
      ADNC.LMH.num.P.Value,
      ADNC.LMH.num.adj.P.Val
    ) %>%
    rename(
      "identifier" = id,
    )

  df1$loc_on_MAPT_8 <- if_else(
    is.na(df1$loc_on_MAPT_8),
    "NA",
    df1$loc_on_MAPT_8
  )
  df1$loc_on_MAPT_8 <- if_else(
    str_detect(df1$loc_on_MAPT_8, "^not"),
    "NA",
    df1$loc_on_MAPT_8
  )
  df1$for_volcano <- if_else(
    df1$loc_on_MAPT_8 == "NA",
    df1$for_volcano,
    df1$loc_on_MAPT_8
  )
  df1 <- df1 %>% select(-loc_on_MAPT_8)
  df1$for_volcano <- str_replace(df1$for_volcano, "Tau441", "TAU")

  df2 <- pannda.sig.ub %>%
    mutate(
      `Gene name` = Gene,
      ADNC.LMH.num.P.Value = `P-Value`,
      ADNC.LMH.num.logFC = -1,
      ADNC.LMH.num.adj.P.Val = FDR,
      identifier = newid,
      site = stringr::str_extract(ModSites, "\\d+"),
      aa = stringr::str_extract(ModSites, "[A-Za-z]+")
    ) %>%
    mutate(
      for_volcano = paste(`Gene name`, site, aa, sep = "_")
    ) %>%
    select(
      identifier,
      for_volcano,
      `Gene name`,
      ADNC.LMH.num.logFC,
      ADNC.LMH.num.P.Value,
      ADNC.LMH.num.adj.P.Val
    ) %>%
    group_by(identifier) %>%
    filter(ADNC.LMH.num.P.Value == min(ADNC.LMH.num.P.Value)) %>%
    ungroup()
  df2$for_volcano <- str_replace(df2$for_volcano, "MAPT", "TAU")
  df <- rbind(df1, df2)
  genes <- df$`Gene name` %>% unique()
}


gene_info <- getBM(
  attributes = c(
    "hgnc_symbol",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand"
  ),
  filters = "hgnc_symbol",
  values = genes,
  mart = ensembl
)

valid_chr <- c(as.character(1:22), "X", "Y", "MT")


gene_info <- gene_info %>%
  dplyr::filter(chromosome_name %in% valid_chr) %>%
  dplyr::select(gene_name = hgnc_symbol, chr = chromosome_name, start_position)


df1 <- df %>% inner_join(gene_info, by = c("Gene name" = "gene_name"))
dim(df1)
df1 <- as.data.frame(df1)


df1$pos <- str_split_i(df1$for_volcano, "_", 2) %>% as.numeric()
df1$chr <- factor(df1$chr, levels = c(as.character(1:22), "X", "Y", "MT"))
df1 <- df1 %>% arrange(chr, start_position, pos)
df1 <- df1 %>%
  group_by(chr) %>%
  mutate(POS = row_number()) %>%
  ungroup()
df1 <- df1 %>% dplyr::select(-pos)
genome_wide_p_val_u <- df1 %>%
  filter(ADNC.LMH.num.logFC > 0, ADNC.LMH.num.adj.P.Val <= 0.05) %>%
  pull(ADNC.LMH.num.P.Value) %>%
  max()
genome_wide_p_val_l <- df1 %>%
  filter(ADNC.LMH.num.logFC < 0, ADNC.LMH.num.adj.P.Val <= 0.05) %>%
  pull(ADNC.LMH.num.P.Value) %>%
  max()


df1 <- df1 %>%
  dplyr::select(
    rsid = for_volcano,
    chr,
    pos = POS,
    beta = ADNC.LMH.num.logFC,
    pval = ADNC.LMH.num.P.Value
  )
df1 <- as.data.frame(df1)

my_upper_colors <- RColorBrewer::brewer.pal(4, "Paired")[2:1]
my_lower_colors <- RColorBrewer::brewer.pal(4, "Paired")[4:3]

plot_df <- df1 %>%
  group_by(chr) %>%
  # Compute chromosome size
  summarise(chrlength = max(pos)) %>%
  # Calculate cumulative position of each chromosome
  mutate(cumulativechrlength = cumsum(as.numeric(chrlength)) - chrlength) %>%
  dplyr::select(-chrlength) %>%
  # Temporarily add the cumulative length of each chromosome to the initial
  # dataset
  left_join(df1, ., by = c("chr" = "chr")) %>%
  # Sort by chr then position
  arrange(chr, pos) %>%
  # Add the position to the cumulative chromosome length to get the position of
  # this probe relative to all other probes
  mutate(rel_pos = pos + cumulativechrlength) %>%
  # Calculate the logged p-value too
  mutate(logged_p = -log10(pval)) %>%
  dplyr::select(-cumulativechrlength)

axis_df <- plot_df %>%
  group_by(chr) %>%
  summarize(
    chr_center = (max(rel_pos) + min(rel_pos)) / 2,
    chr_end = max(rel_pos)
  )
maxp <- ceiling(max(plot_df$logged_p, na.rm = TRUE))


plot_df$colors <- "none"
plot_df$label <- ""


my_upper_colors <- c("#656565", "#bfbfbf")


highlight_size <- 2.3 #
vline_color <- "#6caed5" #
vline_type <- "dotdash" #

plot_df_upper <- plot_df[which(plot_df$beta > 0), ]
plot_df_lower <- plot_df[which(plot_df$beta < 0), ]

up_sig_rsids <- plot_df_upper %>%
  filter(pval <= genome_wide_p_val_u) %>%
  arrange(pval) %>%
  slice_head(n = 30) %>%
  pull(rsid)

lo_sig_rsids <- plot_df_lower %>%
  filter(pval <= genome_wide_p_val_l) %>%
  arrange(pval) %>%
  pull(rsid)


lo_all_rsids <- plot_df_lower %>%
  arrange(pval) %>%
  pull(rsid)

message(
  "top 30 hits in cbmap, detected in pannda ",
  length(unique(intersect(lo_all_rsids, up_sig_rsids)))
)
message(
  "top 30 hits in cbmap, overlapped with pannda top 30 hits ",
  length(unique(intersect(lo_sig_rsids, up_sig_rsids)))
)

highlight_rsids_lw <- intersect(up_sig_rsids, lo_sig_rsids)
highlight_rsids_up <- up_sig_rsids

highlight_data_upper <- plot_df_upper[
  plot_df_upper$rsid %in% highlight_rsids_up,
]
highlight_color_upper <- highlight_data_upper %>%
  mutate(
    pointcolor = if_else(rsid %in% highlight_rsids_lw, "#fff2c2", "#abcfa5")
  ) %>%
  pull(pointcolor)
highlight_data_lower <- plot_df_lower[
  plot_df_lower$rsid %in% highlight_rsids_lw,
]
highlight_color_lower <- highlight_data_lower %>%
  mutate(pointcolor = "#fff2c2") %>%
  pull(pointcolor)

upper_plot <- ggplot(
  data = plot_df_upper,
  aes(x = rel_pos, y = logged_p)
) +
  geom_point(
    aes(color = as.factor(chr)),
    size = 0.8,
    show.legend = F
  ) +
  scale_color_manual(values = rep(my_upper_colors, nrow(axis_df))) +

  geom_segment(
    data = highlight_data_upper,
    aes(
      x = rel_pos, #
      y = logged_p, #
      xend = rel_pos, #
      yend = 0
    ),
    color = vline_color,
    linetype = vline_type,
    linewidth = 0.4
  ) +
  geom_point(
    data = highlight_data_upper,
    aes(x = rel_pos, y = logged_p),
    shape = 21,
    color = "black",
    fill = highlight_color_upper,
    stroke = 0.8,
    size = highlight_size
  ) +
  scale_x_continuous(
    labels = axis_df$chr,
    breaks = axis_df$chr_center,
    expand = expansion(mult = 0.01)
  ) +
  scale_y_continuous(
    limits = c(0, maxp * 1.2),
    expand = expansion(mult = c(0.02, 0))
  ) +
  geom_hline(
    yintercept = -log10(genome_wide_p_val_u),
    color = "red",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  labs(x = "", y = bquote("CBMAP -log"[10] * "(P)")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    plot.margin = margin(t = 10, l = 10, r = 10, b = 0)
  )


lower_plot <- ggplot(
  data = plot_df_lower,
  aes(x = rel_pos, y = logged_p)
) +
  geom_point(
    aes(color = as.factor(chr)),
    size = 0.8,
    show.legend = F
  ) +
  scale_color_manual(values = rep(my_upper_colors, nrow(axis_df))) +

  geom_segment(
    data = highlight_data_lower,
    aes(x = rel_pos, y = logged_p, xend = rel_pos, yend = 0), #
    color = vline_color,
    linetype = vline_type,
    linewidth = 0.4
  ) +

  geom_point(
    data = highlight_data_lower,
    aes(x = rel_pos, y = logged_p),
    shape = 21,

    fill = highlight_color_lower,
    color = "black",
    stroke = 0.8,

    size = highlight_size
  ) +

  scale_x_continuous(
    labels = axis_df$chr,
    breaks = axis_df$chr_center,
    position = "top",
    expand = expansion(mult = 0.01)
  ) +
  scale_y_reverse(
    limits = c(maxp * 1.2, 0),
    expand = expansion(mult = c(0.0, 0.02))
  ) +
  geom_hline(
    yintercept = -log10(genome_wide_p_val_l),
    color = "red",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  labs(x = "", y = bquote("PanNDA -log"[10] * "(P)")) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(t = 0, l = 10, r = 10, b = 10)
  )

p <- upper_plot / lower_plot + plot_layout(guides = "collect")
print(p)


pdf(
  sprintf("../Results/miamiplot/%s_cbmap_vs_pannda.pdf", ptm),
  width = 11,
  height = 5
)
plot(p)
dev.off()


# ==============================================================================
# venn plot
# ==============================================================================

pannda.pr1 <- read_xlsx("../Other_source_data/AD_stats.xlsx", skip = 3) %>%
  mutate(
    ft = str_sub(`Unique identifier`, 1, 2),
    acc = str_split_i(`Unique identifier`, "\\.", 2)
  ) %>%
  filter(
    ft == "sp"
  ) %>%
  rename(id = `Unique identifier`) %>%
  mutate(
    gene_protein = paste0(Gene, "_", acc),
    iso = if_else(
      str_detect(str_split_i(id, "\\.", 3), "^\\d+$"),
      str_split_i(id, "\\.", 3),
      "1"
    )
  ) %>%
  filter(
    iso == "1"
  ) %>%
  select(-iso)
colnames(pannda.pr1) <- c(
  "id",
  "gene",
  "comb_p",
  "fdr",
  "avg_log2FC",
  "avg_log2FC_z",
  "log2FC_1",
  "log2FC_z_1",
  "p_1",
  "log2FC_2",
  "log2FC_z_2",
  "p_2",
  "log2FC_3",
  "log2FC_z_3",
  "p_3",
  "used_in_volcano_plot",
  "ft",
  "acc",
  "gene_protein"
)

pannda.sig.up.pr <- pannda.pr1 %>%
  filter(
    avg_log2FC > 0,
    fdr < 0.05
  ) %>%
  pull(gene) %>%
  unique()

pannda.sig.dn.pr <- pannda.pr1 %>%
  filter(
    avg_log2FC < 0,
    fdr < 0.05
  ) %>%
  pull(gene) %>%
  unique()

pannda.sig.pr <- c(pannda.sig.dn.pr, pannda.sig.up.pr)
length(pannda.sig.pr)

cbmap.sig.up.pr <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 2,
  skip = 1
) %>%
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  pull(genename) %>%
  unique()
cbmap.sig.dn.pr <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 2,
  skip = 1
) %>%
  filter(adj.P.Val < 0.05, logFC < 0) %>%
  pull(genename) %>%
  unique()
cbmap.sig.pr <- c(cbmap.sig.up.pr, cbmap.sig.dn.pr)
length(cbmap.sig.pr)


cbmap.pr <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 2,
  skip = 1
) %>%
  pull(genename) %>%
  unique()
pannda.pr <- read_xlsx("../Other_source_data/AD_stats.xlsx", skip = 3) %>%
  filter(!is.na(Gene)) %>%
  pull(Gene)
setdiff(cbmap.pr, pannda.pr) %>% length()
setdiff(pannda.pr, cbmap.pr) %>% length()

setdiff(cbmap.sig.pr, pannda.pr) %>% length()
setdiff(pannda.sig.pr, cbmap.pr) %>% length()


venn_list <- list(
  "CBMAP sig up" = intersect(cbmap.sig.up.pr, pannda.pr),
  "CBMAP sig dn" = intersect(cbmap.sig.dn.pr, pannda.pr),
  "PanNDA sig up" = intersect(pannda.sig.up.pr, cbmap.pr),
  "PanNDA sig dn" = intersect(pannda.sig.dn.pr, cbmap.pr)
)


v3 <- venn.diagram(
  x = venn_list,
  filename = NULL, #

  fill = brewer.pal(4, "PiYG"), #
  alpha = 0.5, #
  col = brewer.pal(4, "PiYG"), #
  cex = 1.5, #
  cat.cex = 1.2,
  cat.fontface = "bold",
  fontfamily = "Arial",
  fontface = "bold",

  main = ""
)
grid.newpage()
grid.draw(v3)


upset_df <- data.frame(
  gene = unique(c(cbmap.sig.pr, pannda.sig.pr)),
  `Sig up in CBMAP` = F,
  `Sig up in PanNGA` = F,
  `Sig dn in CBMAP` = F,
  `Sig dn in PanNGA` = F,
  `Deteced in another` = "yes",
  check.names = F
)

for (i in seq.int(nrow(upset_df))) {
  tmpgene <- upset_df[i, 1] %>% unlist()
  if (tmpgene %in% cbmap.sig.up.pr) {
    upset_df[i, 2] <- T
  }
  if (tmpgene %in% pannda.sig.up.pr) {
    upset_df[i, 3] <- T
  }
  if (tmpgene %in% cbmap.sig.dn.pr) {
    upset_df[i, 4] <- T
  }
  if (tmpgene %in% pannda.sig.dn.pr) {
    upset_df[i, 5] <- T
  }
  if (
    tmpgene %in%
      c(setdiff(cbmap.sig.pr, pannda.pr), setdiff(pannda.sig.pr, cbmap.pr))
  ) {
    upset_df[i, 6] <- "no"
  }
}
display.brewer.all()
mycols <- brewer.pal(5, "BrBG")
library(ComplexUpset)
genres <- colnames(upset_df)[2:5]
ComplexUpset::upset(
  upset_df,
  genres,
  set_sizes = (upset_set_size() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )),
  base_annotations = list(
    'Intersection size' = intersection_size() +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  ),
  annotations = list(
    'Detected' = (ggplot(mapping = aes(fill = `Deteced in another`)) +
      geom_bar(stat = 'count', position = 'fill') +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_fill_manual(
        values = c(
          "yes" = mycols[1],
          "no" = mycols[2]
        )
      ) +
      ylab('Detected in another dataset') +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ))
  ),
  width_ratio = 0.1
) &
  theme(
    axis.title.x = element_text(
      family = "sans",
      face = "bold",
      color = "transparent",
      size = 0
    ),
    text = element_text(family = "sans", face = "bold")
  )


#===============================================================================
#===============================================================================

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "asia"
) #

if (F) {
  categories <- c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
  gene_sets <- do.call(
    rbind,
    lapply(categories, function(cat) {
      msigdbr(category = cat)
    })
  )
  fwrite(gene_sets, "../Other_source_data/msigdbr_GeneSets_annot.txt")
}

gene_sets <- fread("../Other_source_data/msigdbr_GeneSets_annot.txt")


pannda.pr1 <- read_xlsx("../Other_source_data/AD_stats.xlsx", skip = 3) %>%
  mutate(
    ft = str_sub(`Unique identifier`, 1, 2),
    acc = str_split_i(`Unique identifier`, "\\.", 2)
  ) %>%
  filter(
    ft == "sp"
  ) %>%
  rename(id = `Unique identifier`) %>%
  mutate(
    gene_protein = paste0(Gene, "_", acc),
    iso = if_else(
      str_detect(str_split_i(id, "\\.", 3), "^\\d+$"),
      str_split_i(id, "\\.", 3),
      "1"
    )
  ) %>%
  filter(
    iso == "1"
  ) %>%
  select(-iso)
dim(pannda.pr1)

colnames(pannda.pr1) <- c(
  "id",
  "gene",
  "comb_p",
  "fdr",
  "avg_log2FC",
  "avg_log2FC_z",
  "log2FC_1",
  "log2FC_z_1",
  "p_1",
  "log2FC_2",
  "log2FC_z_2",
  "p_2",
  "log2FC_3",
  "log2FC_z_3",
  "p_3",
  "used_in_volcano_plot",
  "ft",
  "acc",
  "gene_protein"
)


if (F) {
  head(pannda.pr1)

  tmpdata <- pannda.pr1 %>% select(p_1, log2FC_1) %>% filter(p_1 > 0.31)
}


gene_symbols <- unique(pannda.pr1$gene)
mapping <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id", "gene_biotype"),
  filters = "hgnc_symbol",
  values = gene_symbols,
  mart = mart
)
mapping <- mapping[!duplicated(mapping$hgnc_symbol), ]
pannda.pr1$entre_id <- mapping[
  match(pannda.pr1$gene, mapping$hgnc_symbol),
  "entrezgene_id"
]


pannda.pr1.dist <- pannda.pr1 %>%
  filter(!is.na(entre_id), !duplicated(entre_id))
dim(pannda.pr1.dist)


if (F) {
  map <- gene_sets[, c("gs_name", 'entrez_gene')]
  map$entrez_gene <- as.character(map$entrez_gene)
  get_gsea <- function(col) {
    set.seed(2025)
    pannda.pr1.dist <- pannda.pr1.dist %>%
      arrange(desc(.data[[col]]))
    gene_list <- pannda.pr1.dist[[col]]
    names(gene_list) <- pannda.pr1.dist$entre_id
    result_GSEA_lim = GSEA(
      geneList = gene_list,
      TERM2GENE = map,
      eps = 0,
      pvalueCutoff = 1
    ) # 输出所有!
    gseapannda <- result_GSEA_lim@result
    gseapannda <- gseapannda %>% as_tibble()
    gseapannda
  }

  gseapannda0 <- get_gsea(col = "avg_log2FC_z")
  gseapannda1 <- get_gsea(col = "log2FC_z_1")
  gseapannda2 <- get_gsea(col = "log2FC_z_2")
  gseapannda3 <- get_gsea(col = "log2FC_z_3")

  gseapannda0$source <- "PanNDA Avg"
  gseapannda1$source <- "PanNDA Disc"
  gseapannda2$source <- "PanNDA Val1"
  gseapannda3$source <- "PanNDA Val2"
  gseapannda <- Reduce(
    rbind,
    list(gseapannda0, gseapannda1, gseapannda2, gseapannda3)
  )
  fwrite(gseapannda, "../Other_source_data/gsea_pannda_pr_ad.txt")
}


gseapannda <- fread("../Other_source_data/gsea_pannda_pr_ad.txt")
dim(gseapannda)
gseapannda <- gseapannda %>%
  filter(
    source == "PanNDA Avg",
    str_detect(ID, "^(GOBP|GOCC|GOMF|KEGG|REACTOME)_")
  )
nrow(gseapannda)


map.aggt <- gene_sets %>%
  select(gs_name, entrez_gene, ensembl_gene, gene_symbol) %>%
  group_by(gs_name) %>%
  summarise(
    entrez_genes = paste0(entrez_gene, collapse = ","),
    gene_symbols = paste0(gene_symbol, collapse = ",")
  )

gseapannda %>%
  semi_join(
    map.aggt,
    by = c("ID" = "gs_name")
  ) %>%
  nrow()

gseapannda <- gseapannda %>%
  left_join(
    map.aggt,
    by = c("ID" = "gs_name")
  )


dim(gseapannda)


categorize_pathway <- function(description) {
  t <- tolower(description)

  if (
    grepl(
      "splice|snrnp|mrna|ribonucleoprotein|rna processing|spliceosome|transcript|translation|ribosome",
      t
    )
  ) {
    return("RNA Processing & Splicing")
  }

  if (
    grepl(
      "actin|filament|adhesion|matrix|collagen|skeleton|microtubule|myosin|tubulin|focal adhesion|muscle|junction|cytoskeleton|desmosome|binding",
      t
    )
  ) {
    return("Cytoskeleton & Adhesion")
  }

  if (
    grepl(
      "mitochondri|oxidative|respirat|atp|electron|proton|energy|nadh|dehydrogenase|carboxylic acid|tca cycle|fatty acid|metabolic",
      t
    )
  ) {
    return("Mitochondrial & Energy Metabolism")
  }

  if (
    grepl(
      "signaling|pathway|wnt|tgf|mtor|mapk|notch|hedgehog|pi3k|kinase|receptor|cytokine|chemokine|second messenger|g protein|hormonal",
      t
    )
  ) {
    return("Signal Transduction")
  }

  if (
    grepl(
      "differentiation|morphogenesis|development|embryo|osteoblast|myeloid|vasculature|\\borgan\\b|tissue|neurogenesis|chondrocyte",
      t
    )
  ) {
    return("Development & Differentiation")
  }

  if (
    grepl(
      "transport|localization|vesicle|export|import|pore|channel|secretion|trafficking|endocytosis|exocytosis|golgi|\\ber\\b|nucleus|\\bion\\b",
      t
    )
  ) {
    return("Transport & Localization")
  }

  if (
    grepl(
      "response|stress|inflammation|death|apoptotic|starvation|hormone|immune|defense|innate|adaptive|homeostasis|wound",
      t
    )
  ) {
    return("Response & Homeostasis")
  }

  if (
    grepl(
      "disease|cancer|infection|alzheimer|parkinson|huntington|amyloid|malaria|prion",
      t
    )
  ) {
    return("Human Disease")
  }

  return("Others / Cellular Regulation")
}


gseacbmap <- fread("../Other_source_data/gsea_cbmap_pr_adnc.txt")

gseacbmap <- gseacbmap %>%
  filter(
    str_detect(ID, "^(GOBP|GOCC|GOMF|KEGG|REACTOME)_")
  )
gseacbmap
gseacbmap$Category <- sapply(gseacbmap$Description, categorize_pathway)
gseacbmap %>%
  semi_join(
    map.aggt,
    by = c("ID" = "gs_name")
  ) %>%
  nrow()

gseacbmap <- gseacbmap %>%
  left_join(
    map.aggt,
    by = c("ID" = "gs_name")
  )

dim(gseacbmap)
# write.xlsx(gseacbmap, file = "../Results/CBMAP_PanNDA/gseacbmap_with_cat.xlsx")
if (F) {
  goodpathgene.cbmap <- gseacbmap %>%
    filter(Description == "intermediate filament cytoskeleton") %>%
    pull(Gene_Symbols) %>%
    str_split_1(",") %>%
    str_trim()
  goodpathgene.cbmap %>% length()
  goodpathgene.pannda <- gseapannda %>%
    filter(Description == "intermediate filament cytoskeleton") %>%
    pull(gene)
}


gseapannda$ID %>% unique() %>% str_split_i("_", 1) %>% table()
gseacbmap$ID %>% unique() %>% str_split_i("_", 1) %>% table()

gseacbmap <- gseacbmap %>%
  separate(
    col = Description,
    into = c(NA, "Description"),
    sep = "_",
    extra = "merge"
  ) %>%
  mutate(
    Description = str_to_lower(Description) %>% str_replace_all("_", " ")
  )
gseacbmap[, 1:3] %>% head()
gseapannda <- gseapannda %>%
  separate(
    col = Description,
    into = c(NA, "Description"),
    sep = "_",
    extra = "merge"
  ) %>%
  mutate(
    Description = str_to_lower(Description) %>% str_replace_all("_", " ")
  )
gseapannda[, 1:3] %>% head()


gseacbmapsigpath <- gseacbmap %>%
  filter(
    p.adjust < 0.05,
    qvalue < 0.05
  ) %>%
  arrange(Category) %>%
  pull(ID) %>%
  unique()
length(gseacbmapsigpath)
gseacbmappath <- gseacbmap %>% pull(ID) %>% unique()
length(gseacbmappath)

gseapanndasigpath <- gseapannda %>%
  filter(
    p.adjust < 0.05,
    qvalue < 0.05
  ) %>%
  pull(ID) %>%
  unique()
length(gseapanndasigpath)
gseapanndapath <- gseapannda %>% pull(ID) %>% unique()
length(gseapanndapath)


if (F) {
  library(org.Hs.eg.db)
  cbmap.pr <- read_excel(
    "../Other_source_data/result_table.xlsx",
    sheet = 2,
    skip = 1
  ) %>%
    arrange(desc(logFC))

  genesin <- gseapannda %>%
    filter(
      # Description %in%  c("exocytic vesicle", "transport vesicle", "transport vesicle membrane", "synaptic vesicle membrane"),
      str_detect(Description, "spliceosom"), #
      source == "PanNDA Avg"
    ) %>%
    mutate(genes_ls = str_split(gene_symbols, ",")) %>%
    select(genes_ls) %>%
    unnest(cols = c(genes_ls)) %>%
    distinct(genes_ls) %>%
    left_join(
      pannda.pr1[, c("gene", "log2FC_1", "comb_p", "fdr")],
      by = c("genes_ls" = "gene")
    ) %>%
    left_join(
      cbmap.pr[, c("logFC", "P.Value", "adj.P.Val", "genename")],
      by = c("genes_ls" = "genename")
    ) %>%
    filter(
      sign(log2FC_1) != sign(logFC)
    ) %>%
    arrange(
      adj.P.Val
    )
  genesin %>% print(n = Inf)

  venn_list <- list(
    "CBMAP path" = gseacbmappath,
    "PanNDA path" = gseapanndapath,
    "CBMAP sig path" = gseacbmapsigpath,
    "PanNDA sig path" = gseapanndasigpath
  )

  library(grid)

  v3 <- venn.diagram(
    x = venn_list,
    filename = NULL, #

    fill = brewer.pal(4, "PiYG"),
    alpha = 0.5, #
    col = brewer.pal(4, "PiYG"), #

    cex = 1.5, #
    cat.cex = 1.5, #
    cat.fontface = "bold", #
    fontfamily = "Arial",
    fontface = "bold",

    main = ""
  )

  grid.newpage()
  grid.draw(v3)
}


gseacbmapsig <- gseacbmap %>%
  filter(
    p.adjust < 0.05,
    qvalue < 0.05
  )
dim(gseacbmapsig)
# write.xlsx(gseacbmapsig, file = "../Results/CBMAP_PanNDA/gseacbmapsig.xlsx")

gseapannda.bak <- gseapannda
gseapannda <- gseapannda.bak
gseapannda <- gseapannda %>%
  right_join(
    gseacbmapsig[, c("ID", "Description", "Category")],
    by = c("ID" = "ID", "Description" = "Description")
  )

dim(gseapannda)
head(gseapannda[, 1:10], 10)
sum(is.na(gseapannda$Description))
sum(is.na(gseapannda$NES))

gseapannda <- gseapannda %>%
  mutate(
    NES = if_else(is.na(NES), 0, NES),
    # pvalue = if_else(is.na(pvalue), 0, pvalue),
    # qvalue = if_else(is.na(qvalue), 0, qvalue),
    # p.adjust = if_else(is.na(p.adjust), 0, p.adjust),
    source = if_else(is.na(source), "PanNDA Avg", source)
  )


# gseacbmapsig %>% select(Description, Category) %>% arrange(Category) %>% write.xlsx("tmpcbmap.xlsx")

head(gseacbmapsig, 10) %>% as_tibble()
cat_genes <- gseacbmapsig %>%
  group_by(Category) %>%
  summarise(
    gene = paste(na.omit(gene_symbols), collapse = ",")
  )

cat2gsea <- data.frame()

for (i in seq.int(nrow(cat_genes))) {
  gns <- cat_genes[i, ]$gene
  gns <- str_split_1(gns, ",|_") %>% str_trim(side = "both") %>% unique()
  tmp <- data.frame(
    term = cat_genes[i, ]$Category,
    gene = gns
  )
  cat2gsea <- rbind(cat2gsea, tmp)
}
cat2gsea <- cat2gsea %>% distinct()
cat2gsea %>% head()
dim(cat2gsea)


cats <- c(
  "Others / Cellular Regulation",
  "Mitochondrial & Energy Metabolism",
  "Human Disease",
  "Transport & Localization",
  "RNA Processing & Splicing",
  "Signal Transduction",
  "Cytoskeleton & Adhesion",
  "Response & Homeostasis",
  "Development & Differentiation"
)


gseacbmapsig.dedup <- gseacbmapsig %>%
  mutate(source = "CBMAP") %>%
  group_by(Description) %>%
  slice_head(n = 1) %>%
  ungroup()
gseapannda.dedup <- gseapannda %>%
  group_by(Description) %>%
  slice_head(n = 1) %>%
  ungroup()

gseacbmapsig.ht <- gseacbmapsig.dedup %>%
  select(
    Description,
    Category,
    source,
    NES,
    p.adjust
  )
gseapannda.ht <- gseapannda.dedup %>%
  select(
    Description,
    Category,
    source,
    NES,
    p.adjust
  )


gseapath <- rbind(gseacbmapsig.ht, gseapannda.ht) %>%
  pivot_wider(
    id_cols = Description:Category,
    names_from = source,
    values_from = NES
  ) %>%
  filter(
    Category %in% cats
  ) %>%
  mutate(
    Category = factor(Category, levels = cats)
  )

gseapval <- rbind(gseacbmapsig.ht, gseapannda.ht) %>%
  pivot_wider(
    id_cols = Description:Category,
    names_from = source,
    values_from = p.adjust
  ) %>%
  filter(
    Category %in% cats
  ) %>%
  mutate(
    Category = factor(Category, levels = cats)
  )


tmp.res <- data.frame()

for (i in cats) {
  tmp <- gseapath %>%
    filter(Category == i)
  cluster_data <- tmp[, seq.int(3, ncol(gseapath))]
  cluster_data_scaled <- scale(cluster_data) #
  d <- dist(cluster_data_scaled)
  hc <- hclust(d, method = "complete")
  tmp <- tmp[hc$order, ]
  tmp.res <- rbind(tmp.res, tmp)
}
tmp.res
gseapath <- tmp.res %>%
  column_to_rownames(
    "Description"
  ) %>%
  select(
    -ends_with("y")
  )
gseapval <- gseapval[match(rownames(gseapath), gseapval$Description), ] %>%
  column_to_rownames(
    "Description"
  ) %>%
  select(
    -ends_with("y")
  )
split <- tmp.res$Category


library(RColorBrewer)

display.brewer.all()
mycolors <- brewer.pal(9, "Set3")

gseapath <- as.matrix(gseapath)
tmp <- gseapath
rownames(tmp) <- NULL

# sig_color <- gseapannda.dedup %>%
#   filter(
#     source == "PanNDA Avg",
#     Description %in% rownames(gseapath)) %>%
#   mutate(
#     sig = if_else(qvalue < 0.05 & p.adjust < 0.05, "red", "black")
#   ) %>%
#   select(
#     Description, sig
#   ) %>%
#   as.data.frame() %>%
#   column_to_rownames("Description") %>%
#   as.matrix()
#
# sig_color <- sig_color[, 1]
# sig_color <- sig_color[rownames(gseapath)]
ha_row <- rowAnnotation(
  Category = split,
  col = list(
    Category = c(
      "RNA Processing & Splicing" = "#C3D9F0",
      "Cytoskeleton & Adhesion" = "#FFF2C2",
      "Mitochondrial & Energy Metabolism" = "#B5E6DE",
      "Transport & Localization" = "#D3C9EB",
      "Response & Homeostasis" = "#F9D8D6",
      "Development & Differentiation" = mycolors[6],
      "Human Disease" = mycolors[7],
      "Signal Transduction" = mycolors[8],
      "Others / Cellular Regulation" = mycolors[9]
    )
  ),
  Genes = anno_text(
    rownames(gseapath),
    gp = gpar(
      col = "black",
      fontfamily = "sans",
      fontface = "italic",
      fontsize = 10 #
    ),
    just = "left",
    location = unit(0, "npc")
  ),
  show_annotation_name = F,
  annotation_name_gp = gpar(fontsize = 16),
  # annotation_width = unit(rep(8, 6), "mm"),
  simple_anno_size_adjust = T,
  # width = unit(39, "mm"),

  # simple_anno_size = unit(15, "mm"),
  annotation_name_rot = 45,
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 9)
  )
)

mycols <- brewer.pal(8, "Accent")
ha_col <- HeatmapAnnotation(
  # Source = c("CBMAP", "PanNDA", "PanNDA", "PanNDA", "PanNDA"),
  Source = c("CBMAP", "PanNDA"),
  col = list(
    PTM = c("CBMAP" = mycols[5], "PanNDA" = mycols[6])
  ),
  show_annotation_name = F,
  annotation_name_gp = gpar(fontsize = 16),
  # annotation_height = unit(rep(8, 2), "mm"),
  simple_anno_size_adjust = TRUE,
  gp = gpar(col = "darkgrey", lwd = .5),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 9)
  )
)


ht <- Heatmap(
  tmp,
  row_split = split,
  row_title = NULL,
  name = "CBMAP PanNDA",
  col = circlize::colorRamp2(
    c(min(gseapath), 0, max(gseapath)),
    c("#3C5488FF", "white", "#DC0000FF")
  ),
  cluster_rows = T,
  cluster_columns = FALSE, #

  show_row_dend = T,
  show_column_dend = FALSE,

  show_column_names = T,
  show_row_names = T,

  right_annotation = ha_row,
  bottom_annotation = ha_col,
  # heatmap_height = unit(17, "inches"),
  # heatmap_width = unit(10, "inches"),   #
  # row_dend_width = unit(1, "inches"),
  heatmap_legend_param = list(
    title = "Normalized Enrichment Score",
    title_gp = gpar(fontsize = 12),
    labels_gp = gpar(fontsize = 9),
    direction = "horizontal"
  ),
  column_names_rot = 45,
  width = unit(4, "inches"), #
  height = unit(28, "inches"),
  border = FALSE
)


pdf("../Results/CBMAP_PanNDA/pathway_htmp_1225.pdf", 14, 35)
# png("../Results/CBMAP_PanNDA/pathway_htmp_1225.png", 14, 23, units = "in", res = 900)
draw(
  ht,
  padding = unit(c(2, -40, 2, 20), "mm"),
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE
)
dev.off()


for (i in cats) {
  subtmp <- tmp[which(split == i), ]
  subgseapath <- gseapath[which(split == i), ]
  subgseapval <- gseapval[which(split == i), ]
  # subgseapval[, 1] <- 1
  ha_row <- rowAnnotation(
    Category = split[split == i],
    col = list(
      Category = c(
        "RNA Processing & Splicing" = "#C3D9F0",
        "Cytoskeleton & Adhesion" = "#FFF2C2",
        "Mitochondrial & Energy Metabolism" = "#B5E6DE",
        "Transport & Localization" = "#D3C9EB",
        "Response & Homeostasis" = "#F9D8D6",
        "Development & Differentiation" = mycolors[6],
        "Human Disease" = mycolors[7],
        "Signal Transduction" = mycolors[8],
        "Others / Cellular Regulation" = mycolors[9]
      )
    ),
    Genes = anno_text(
      rownames(subgseapath),
      gp = gpar(
        col = "black",
        fontfamily = "sans",
        fontface = "italic",
        fontsize = 10 #
      ),
      just = "left",
      location = unit(0, "npc")
    ),
    show_annotation_name = F,
    annotation_name_gp = gpar(fontsize = 16),
    # annotation_width = unit(rep(8, 6), "mm"),
    simple_anno_size_adjust = T,
    # width = unit(39, "mm"),

    # simple_anno_size = unit(15, "mm"),
    annotation_name_rot = 45,
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 9)
    )
  )
  colfn <- circlize::colorRamp2(
    c(min(subtmp), 0, max(subtmp)),
    c("#3C5488FF", "white", "#DC0000FF")
  )
  if (min(subtmp) >= 0) {
    colfn <- circlize::colorRamp2(
      c(min(subtmp), max(subtmp)),
      c("white", "#DC0000FF")
    )
  }
  if (max(subtmp) <= 0) {
    colfn <- circlize::colorRamp2(
      c(min(subtmp), max(subtmp)),
      c("#3C5488FF", "white")
    )
  }
  ht <- Heatmap(
    subtmp,
    row_title = NULL,
    name = "CBMAP PanNDA",
    col = colfn,
    cluster_rows = T,
    cluster_columns = FALSE, #

    show_row_dend = T,
    show_column_dend = FALSE,

    show_column_names = T,
    # show_row_names = T,

    right_annotation = ha_row,
    bottom_annotation = ha_col,
    # row_dend_width = unit(1, "inches"),
    heatmap_legend_param = list(
      title = "Normalized Enrichment Score",
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 9),
      direction = "horizontal"
    ),
    layer_fun = function(j, i, x, y, width, height, fill) {
      p_values <- subgseapval[cbind(i, j)]

      stars <- sapply(p_values, function(p) {
        if (is.na(p)) {
          return("NA")
        }
        if (p < 0.00000001) {
          return("********")
        }
        if (p < 0.0000001) {
          return("*******")
        }
        if (p < 0.000001) {
          return("******")
        }
        if (p < 0.00001) {
          return("*****")
        }
        if (p < 0.0001) {
          return("****")
        }
        if (p < 0.001) {
          return("***")
        }
        if (p < 0.01) {
          return("**")
        }
        if (p < 0.05) {
          return("*")
        }
        return("")
      })
      vjust_vec <- if_else(stars == "NA", 0.5, 0.75)
      grid.text(
        stars,
        x,
        y,
        vjust = vjust_vec,
        gp = gpar(
          fontsize = 12,
          col = "#2C3E50", #
          fontface = "bold"
        )
      )
    },
    column_names_rot = 45,
    width = unit(3, "inches"), #
    height = unit(15, "inches"),
    border = FALSE
  )
  pdf(
    sprintf(
      "../Results/CBMAP_PanNDA/pathway_sub_htmp_%s_1225.pdf",
      str_replace_all(i, "/", "&")
    ),
    18,
    18
  )
  # png(sprintf("../Results/CBMAP_PanNDA/pathway_sub_htmp_%s_1225.png", i), 8, 12, units = "in", res = 900)
  draw(
    ht,
    padding = unit(c(2, 20, 2, 20), "mm"),
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE
  )
  dev.off()
}


#===============================================================================
# heatmap showing pval between PanNDA and CBMAP
#===============================================================================
cbmap.sig.pr <- read_excel(
  "../Other_source_data/result_table.xlsx",
  sheet = 2,
  skip = 1
) %>%
  arrange(desc(logFC)) %>%
  filter(
    adj.P.Val < 0.05
  )

cbmap.sig.pr.cat <- read_xlsx(
  "../Other_source_data/CBMAP_sigPrs_category_tuning.xlsx",
  sheet = 2
)

if (F) {
  # fine-tuning gene category.......
  table(cbmap.sig.pr.cat$Category)

  cbmap.sig.pr.cat <- cbmap.sig.pr.cat %>%
    mutate(
      Category = case_when(
        Category == "Cytoskeleton & Adhesion" ~ "Cytoskeleton & Myelin",
        Category == "Mitochondrial & Energy Metabolism" ~ "Energy",
        Category == "RNA Processing & Splicing" ~ "RNA splicing",
        Category == "Transport & Localization" ~ "Synaptic vesicle cycle",
        T ~ Category
      )
    ) %>%
    mutate(
      Category = if_else(
        gene %in%
          c(
            "CTHRC1",
            "MFGE8",
            "COL6A2",
            "EFEMP2",
            "FBLN5",
            "EDIL3",
            "DCN",
            "COL6A3"
          ),
        "Extracellular matrix",
        Category
      )
    )
  table(cbmap.sig.pr.cat$Category)
  write_xlsx(cbmap.sig.pr.cat, "../Other_source_data/tmp.xlsx")
}


cbmap.sig.pr <- cbmap.sig.pr %>%
  inner_join(cbmap.sig.pr.cat, by = c("genename" = "gene"))
n_distinct(cbmap.sig.pr$genename) #205

pannda.pr1 <- read_xlsx("../Other_source_data/AD_stats.xlsx", skip = 3) %>%
  mutate(
    ft = str_sub(`Unique identifier`, 1, 2),
    acc = str_split_i(`Unique identifier`, "\\.", 2)
  ) %>%
  filter(
    ft == "sp"
  ) %>%
  rename(id = `Unique identifier`) %>%
  mutate(
    gene_protein = paste0(Gene, "_", acc),
    iso = if_else(
      str_detect(str_split_i(id, "\\.", 3), "^\\d+$"),
      str_split_i(id, "\\.", 3),
      "1"
    )
  ) %>%
  filter(
    iso == "1"
  ) %>%
  select(-iso)


colnames(pannda.pr1) <- c(
  "id",
  "gene",
  "comb_p",
  "fdr",
  "avg_log2FC",
  "avg_log2FC_z",
  "log2FC_1",
  "log2FC_z_1",
  "p_1",
  "log2FC_2",
  "log2FC_z_2",
  "p_2",
  "log2FC_3",
  "log2FC_z_3",
  "p_3",
  "ft",
  "acc",
  "gene_protein"
)


pannda.pr1 <- pannda.pr1 %>%
  filter(gene %in% cbmap.sig.pr$genename) %>%
  inner_join(cbmap.sig.pr.cat, by = "gene") %>%
  distinct(gene, .keep_all = TRUE)

dim(pannda.pr1)


# insoluble proteome

pannda.pr2 <- read_csv("../Other_source_data/AD_insoluble.csv") %>%
  mutate(
    ft = str_sub(Identifier, 1, 2),
    acc = str_split_i(Identifier, "\\.", 2)
  ) %>%
  filter(
    ft == "sp"
  ) %>%
  rename(id = Identifier) %>%
  mutate(
    gene_protein = paste0(Gene, "_", acc),
    iso = if_else(
      str_detect(str_split_i(id, "\\.", 3), "^\\d+$"),
      str_split_i(id, "\\.", 3),
      "1"
    )
  ) %>%
  filter(
    iso == "1"
  ) %>%
  select(-iso)
colnames(pannda.pr2) <- c(
  "id",
  "gene",
  "enrich_factor",
  "p",
  "fdr",
  "log2FC",
  "log2FC_z",
  "ft",
  "acc",
  "gene_protein"
)

pannda.pr2 <- pannda.pr2 %>%
  filter(gene %in% cbmap.sig.pr$genename) %>%
  inner_join(cbmap.sig.pr.cat, by = "gene") %>%
  distinct(gene, .keep_all = T)

dim(pannda.pr2)


pannda.pr1 <- pannda.pr1 %>%
  select(gene, Category, fdr, comb_p) %>%
  mutate(source = "PanNDA Whole") %>%
  rename(new_p = comb_p)

pannda.pr2 <- pannda.pr2 %>%
  select(gene, Category, fdr, p) %>%
  mutate(source = "PanNDA Inslb")

cbmap.sig.pr <- cbmap.sig.pr %>%
  select(genename, P.Value, Category, adj.P.Val) %>%
  mutate(source = "CBMAP")


# sub.cats <- c(
#   "RNA Processing & Splicing",
#   "Cytoskeleton & Adhesion",
#   "Mitochondrial & Energy Metabolism",
#   "Transport & Localization",
#   "Response & Homeostasis"
# )
sub.cats <- c(
  "RNA splicing",
  "Cytoskeleton & Myelin",
  "Energy",
  "Extracellular matrix",
  "Synaptic vesicle cycle"
)


pval.df <- cbmap.sig.pr %>%
  inner_join(
    pannda.pr1,
    by = c("genename" = "gene", "Category" = "Category")
  ) %>%
  left_join(
    pannda.pr2,
    by = c("genename" = "gene", "Category" = "Category")
  ) %>%
  select(-starts_with("source")) %>%
  rename(
    "CBMAP" = "P.Value",
    "PanNDA Inslb" = "p",
    "PanNDA whole" = "new_p"
  ) %>%
  filter(
    Category %in% sub.cats
  ) %>%
  arrange(
    Category
  )
dim(pval.df)
known.AD.pr <- read_csv(
  "../Other_source_data/alzheimer_gene_literature_final.csv"
)
known.AD.pr <- known.AD.pr %>%
  filter(
    gene_name %in% pval.df$genename,
    type %in% c("association and function", "function"),
    str_detect(key_word, "animal")
  )
known.AD.pr


split <- pval.df$Category
pval.mat <- pval.df %>%
  select(-Category, -starts_with("fdr"), -adj.P.Val) %>%
  column_to_rownames("genename") %>%
  as.matrix()
pval.mat <- -log10(pval.mat)
dim(pval.mat)
fdr.mat <- pval.df %>%
  select(-Category, -CBMAP, -starts_with("PanNDA")) %>%
  column_to_rownames("genename") %>%
  rename(
    "CBMAP" = "adj.P.Val",
    "PanNDA whole" = "fdr.x",
    "PanNDA Inslb" = "fdr.y"
  ) %>%
  as.matrix()


mycolors <- brewer.pal(9, "Set3")


mycols <- brewer.pal(8, "Accent")


if (F) {
  ha_row <- rowAnnotation(
    Category = split,
    col = list(
      Category = c(
        "RNA splicing" = "#C3D9F0",
        "Cytoskeleton & Myelin" = "#FFF2C2",
        "Energy" = "#B5E6DE",
        "Synaptic vesicle cycle" = "#D3C9EB",
        "Extracellular matrix" = "#F9D8D6",
        "Development & Differentiation" = mycolors[6],
        "Human Disease" = mycolors[7],
        "Signal Transduction" = mycolors[8],
        "Others / Cellular Regulation" = mycolors[9]
      )
    ),
    Genes = anno_text(
      rownames(pval.mat),
      gp = gpar(
        col = "black",
        fontfamily = "sans",
        fontface = "italic",
        fontsize = 10 # 设置字号
      ),
      just = "left",
      location = unit(0, "npc")
    ),
    show_annotation_name = F,
    annotation_name_gp = gpar(fontsize = 16),
    # annotation_width = unit(rep(8, 6), "mm"),
    simple_anno_size_adjust = T,
    # width = unit(39, "mm"),
    # simple_anno_size = unit(15, "mm"),
    annotation_name_rot = 45,
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 9)
    )
  )

  mycols <- brewer.pal(8, "Accent")
  ha_col <- HeatmapAnnotation(
    # Source = c("CBMAP", "PanNDA", "PanNDA", "PanNDA", "PanNDA"),
    Source = c("CBMAP", "PanNDA whole", "PanNDA Inslb"),
    col = list(
      PTM = c(
        "CBMAP" = mycols[5],
        "PanNDA whole" = mycols[6],
        "PanNDA Inslb" = mycols[7]
      )
    ),
    show_annotation_name = F,
    annotation_name_gp = gpar(fontsize = 16),
    # annotation_height = unit(rep(8, 2), "mm"),
    simple_anno_size_adjust = TRUE,
    gp = gpar(col = "darkgrey", lwd = .5),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 9)
    )
  )

  ht <- Heatmap(
    pval.mat,
    row_split = split,
    row_title = NULL,
    name = "CBMAP PanNDA",
    col = circlize::colorRamp2(
      c(min(pval.mat, na.rm = T), 0, max(pval.mat, na.rm = T)),
      c("#3C5488FF", "white", "#DC0000FF")
    ),
    cluster_rows = T,
    cluster_columns = FALSE, # 都不聚类

    show_row_dend = T,
    show_column_dend = FALSE,

    show_column_names = T,
    show_row_names = F,

    right_annotation = ha_row,
    bottom_annotation = ha_col,
    # heatmap_height = unit(17, "inches"),
    # heatmap_width = unit(10, "inches"),   #
    # row_dend_width = unit(1, "inches"),
    heatmap_legend_param = list(
      title = "Normalized Enrichment Score",
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 9),
      direction = "horizontal"
    ),
    column_names_rot = 45,
    width = unit(4, "inches"), #
    height = unit(28, "inches"),
    border = FALSE
  )

  pdf("../Results/CBMAP_PanNDA/pval_htmp_0104.pdf", 14, 35)
  draw(
    ht,
    padding = unit(c(2, -40, 2, 20), "mm"),
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE
  )
  dev.off()
}


for (i in sub.cats) {
  subtmp <- pval.mat[which(split == i), ]

  subfdr <- fdr.mat[which(split == i), ]

  ha_row <- rowAnnotation(
    Category = split[split == i],
    col = list(
      Category = c(
        "RNA splicing" = "#C3D9F0",
        "Cytoskeleton & Myelin" = "#FFF2C2",
        "Energy" = "#B5E6DE",
        "Synaptic vesicle cycle" = "#D3C9EB",
        "Extracellular matrix" = "#F9D8D6",
        "Development & Differentiation" = mycolors[6],
        "Human Disease" = mycolors[7],
        "Signal Transduction" = mycolors[8],
        "Others / Cellular Regulation" = mycolors[9]
      )
    ),
    Genes = anno_text(
      rownames(subtmp),
      gp = gpar(
        col = "black",
        fontfamily = "sans",
        fontface = "italic",
        fontsize = 10 #
      ),
      just = "left",
      location = unit(0, "npc")
    ),
    show_legend = F,
    show_annotation_name = F,
    annotation_name_gp = gpar(fontsize = 16),
    # annotation_width = unit(rep(8, 6), "mm"),
    simple_anno_size_adjust = T,
    # width = unit(39, "mm"),

    # simple_anno_size = unit(15, "mm"),
    annotation_name_rot = 45,
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 9)
    )
  )

  colfn <- circlize::colorRamp2(
    c(0, max(subtmp, na.rm = T)),
    c("white", "#DC0000FF")
  )

  ht <- Heatmap(
    subtmp,
    row_title = NULL,
    name = "CBMAP PanNDA",
    col = colfn,
    cluster_rows = T,
    cluster_columns = FALSE, #

    show_row_dend = T,
    show_column_dend = FALSE,

    show_column_names = T,
    show_row_names = F,

    right_annotation = ha_row,
    heatmap_legend_param = list(
      title = bquote("-log"[10] * "(P)"),
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 9),
      direction = "vertical"
    ),
    layer_fun = function(j, i, x, y, width, height, fill) {
      fdrs <- subfdr[cbind(i, j)]
      ps <- subtmp[cbind(i, j)]
      stars <- mapply(
        function(fdr, p) {
          if (is.na(p)) {
            return("NA")
          }
          if (fdr < 0.05) {
            return("**")
          }
          if (p > -log10(0.05)) {
            return("*")
          }
          return("")
        },
        fdrs,
        ps
      )
      vjust_vec <- if_else(stars == "NA", 0.5, 0.75)
      grid.text(
        stars,
        x,
        y,
        vjust = vjust_vec,
        gp = gpar(
          fontsize = 8,
          col = "#2C3E50", #
          fontface = "bold"
        )
      )
    },
    column_names_rot = 45,
    width = unit(1.2, "inches"), #
    height = unit(0.18 * nrow(subtmp), "inches"),
    border = FALSE
  )
  pdf(
    sprintf(
      "../Results/CBMAP_PanNDA/pval_sub_htmp_%s_0114.pdf",
      str_replace_all(i, "/", "&")
    ),
    width = 5,
    height = 3 + 0.18 * nrow(subtmp)
  )
  # png(sprintf("../Results/CBMAP_PanNDA/pathway_sub_htmp_%s_1225.png", i), 8, 12, units = "in", res = 900)
  draw(
    ht,
    padding = unit(c(2, 2, 2, 2), "mm"),
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legend = TRUE
  )
  dev.off()
}
