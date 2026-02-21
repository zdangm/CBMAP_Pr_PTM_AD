suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readxl)
  library(data.table)
  library(biomaRt)
  library(ggplot2)
  library(ggrepel)
  library(gridExtra)
  library(grid)
  library(ACAT)
  library(writexl)
})
rm(list = ls())
conflicted::conflict_prefer_all("dplyr")
setwd("/share/home/lik/Scripts/likun/Random_tasks/Results/miamiplot/")
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "www"
)


files <- list(
  phospho = "/data/projects/China_Brain_MultiOmics/humanBrain_Phospho/PhosPho_maintained_by_LK/DEA/Results/version_after_2025_05_23/log_impute_scale_INT/All.sites.All.traits.intensity.limma.xlsx",
  ubiq = "/data/projects/China_Brain_MultiOmics/humanBrain_Ubiquitylation/Ubiquity_maintained_by_LK/DEA/Results/version_after_2025_05_23/log_impute_scale_INT/All.sites.All.traits.intensity.limma.xlsx",
  ace = "/data/projects/China_Brain_MultiOmics/humanBrain_Acetylation/Acetylation_maintained_by_LK/DEA/Results/version_after_2025_06_30/log_impute_scale_INT/All.sites.All.traits.intensity.limma.xlsx"
)


ptm <- "phospho"
ptm <- "ubiq"
ptm <- "ace"


df <- read_xlsx(files[[ptm]])
colnames(df)[1] <- "identifier"

df <- df %>%
  dplyr::select(
    identifier,
    `Gene name`,
    loc_on_MAPT_8,
    for_volcano,
    braak.num.logFC,
    braak.num.P.Value,
    braak.num.adj.P.Val,
    ADNC.LMH.num.logFC,
    ADNC.LMH.num.P.Value,
    ADNC.LMH.num.adj.P.Val,
    adj3.ADNC.LMH.num.logFC,
    adj3.ADNC.LMH.num.P.Value,
    adj3.ADNC.LMH.num.adj.P.Val,
    a.score.logFC,
    a.score.P.Value,
    a.score.adj.P.Val,
    b.score.logFC,
    b.score.P.Value,
    b.score.adj.P.Val
  )
df.bak <- df


#-------------------------------here
df <- df.bak
must_cols <- c("Gene name")
ad.trait <- "lmh"

if (ad.trait == "b6") {
  df <- df %>%
    dplyr::select(all_of(must_cols), starts_with("braak")) %>%
    rename(
      ADNC.LMH.num.logFC = braak.num.logFC,
      ADNC.LMH.num.P.Value = braak.num.P.Value,
      ADNC.LMH.num.adj.P.Val = braak.num.adj.P.Val
    )
} else if (ad.trait == "lmh") {
  df <- df %>% dplyr::select(all_of(must_cols), starts_with("ADNC.LMH"))
} else if (ad.trait == "adj.lmh") {
  df <- df %>%
    dplyr::select(all_of(must_cols), starts_with("adj.ADNC.LMH")) %>%
    rename(
      ADNC.LMH.num.logFC = adj.ADNC.LMH.num.logFC,
      ADNC.LMH.num.P.Value = adj.ADNC.LMH.num.P.Value,
      ADNC.LMH.num.adj.P.Val = adj.ADNC.LMH.num.adj.P.Val
    )
} else if (ad.trait == "a3") {
  df <- df %>%
    dplyr::select(all_of(must_cols), starts_with("a.score")) %>%
    rename(
      ADNC.LMH.num.logFC = a.score.logFC,
      ADNC.LMH.num.P.Value = a.score.P.Value,
      ADNC.LMH.num.adj.P.Val = a.score.adj.P.Val
    )
} else if (ad.trait == "b3") {
  df <- df %>%
    dplyr::select(all_of(must_cols), starts_with("b.score")) %>%
    rename(
      ADNC.LMH.num.logFC = b.score.logFC,
      ADNC.LMH.num.P.Value = b.score.P.Value,
      ADNC.LMH.num.adj.P.Val = b.score.adj.P.Val
    )
}

colnames(df)
dim(df)

{
  df.acat <- df %>%
    group_by(`Gene name`) %>%
    summarise(
      acat_p_adnc_lmh = ACAT(ADNC.LMH.num.P.Value),
      same_dir = all(ADNC.LMH.num.logFC > 0) | all(ADNC.LMH.num.logFC <= 0)
    ) %>%
    arrange(acat_p_adnc_lmh)
  dim(df.acat)
  df.acat$adj_acat_p_adnc_lmh <- p.adjust(
    df.acat$acat_p_adnc_lmh,
    method = "BH"
  )
  # df.acat %>% write_xlsx(path = sprintf("ACAT_res_%s.xlsx", ptm))

  ori.genes <- df %>%
    filter(ADNC.LMH.num.adj.P.Val < 0.05) %>%
    pull(`Gene name`) %>%
    unique()
  acat.genes <- df.acat %>%
    filter(adj_acat_p_adnc_lmh < 0.05) %>%
    pull(`Gene name`) %>%
    unique()
  length(ori.genes)
  length(acat.genes)
  intersect(ori.genes, acat.genes) %>% length()
  setdiff(acat.genes, ori.genes) %>% length()
  paste(acat.genes, collapse = "','")

  genes <- df.acat$`Gene name` %>% unique()
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
    dplyr::select(
      gene_name = hgnc_symbol,
      chr = chromosome_name,
      start_position
    )

  df1 <- df.acat %>% inner_join(gene_info, by = c("Gene name" = "gene_name"))
  dim(df1)
  df1 <- as.data.frame(df1)
  head(df1)
  df1$chr <- factor(df1$chr, levels = c(as.character(1:22), "X", "Y", "MT"))

  df1 <- df1 %>% arrange(chr, start_position)
  df1 <- df1 %>%
    group_by(chr) %>%
    mutate(POS = row_number()) %>%
    ungroup()

  # df1$acat_p_adnc_lmh <- -log10(df1$acat_p_adnc_lmh)
  genome_wide_p_val <- df1 %>%
    filter(adj_acat_p_adnc_lmh <= 0.05) %>%
    pull(acat_p_adnc_lmh) %>%
    max()

  df1 <- df1 %>%
    dplyr::select(
      rsid = `Gene name`,
      chr,
      pos = POS,
      beta = 1,
      pval = acat_p_adnc_lmh
    )
  df1 <- as.data.frame(df1)

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
    mutate(logged_p = -log10(pval)) %>% # D.Z recommends log10log10, 不现实
    dplyr::select(-cumulativechrlength)

  axis_df <- plot_df %>%
    group_by(chr) %>%
    summarize(
      chr_center = (max(rel_pos) + min(rel_pos)) / 2,
      chr_end = max(rel_pos)
    )
  maxp <- ceiling(max(plot_df$logged_p, na.rm = TRUE))

  useTOP <- T

  if (useTOP) {
    message("use top hits")

    top_hits <- plot_df %>%
      filter(beta > 0, pval <= genome_wide_p_val) %>%
      arrange(pval) %>%
      slice_head(n = 10)

    message(nrow(top_hits), " top hits")
  }

  fine_tuning.cat <- NULL

  if (ptm == "phospho") {
    fine_tuning.cat <- read_excel(
      "/data/projects/lik/CBMAP/Random_tasks/Other_source_data/CBMAP_sigPrs_category_tuning.xlsx",
      sheet = "ph"
    )
  } else if (ptm == "ubiq") {
    fine_tuning.cat <- read_excel(
      "/data/projects/lik/CBMAP/Random_tasks/Other_source_data/CBMAP_sigPrs_category_tuning.xlsx",
      sheet = "ub"
    )
  } else if (ptm == "ace") {
    fine_tuning.cat <- read_excel(
      "/data/projects/lik/CBMAP/Random_tasks/Other_source_data/CBMAP_sigPrs_category_tuning.xlsx",
      sheet = "ac"
    )
  }

  plot_df <- plot_df %>%
    left_join(
      fine_tuning.cat,
      by = c("rsid" = "gene")
    ) %>%
    mutate(
      Category = if_else(is.na(Category), "Other", Category)
    ) %>%
    rename(colors = Category)

  col_value <- c(
    "RNA splicing" = "#C3D9F0",
    "Cytoskeleton & Myelin" = "#FFF2C2",
    "Energy metabolism" = "#B5E6DE",
    "Synaptic vesicle cycle" = "#D3C9EB",
    "Extracellular matrix" = "#F9D8D6",
    "Hormonal and circadian regualtion" = "#E8DCCF",
    "Other" = "white"
  )

  if (ptm == "phospho") {
    col_value <- c(
      "RNA splicing" = "#C3D9F0",
      "Cytoskeleton & Myelin" = "#FFF2C2",
      "Energy metabolism" = "white",
      "Synaptic vesicle cycle" = "white",
      "Extracellular matrix" = "white",
      "Hormonal and circadian regualtion" = "white",
      "Other" = "white"
    )
  } else if (ptm == "ubiq") {
    col_value <- c(
      "RNA splicing" = "white",
      "Cytoskeleton & Myelin" = "#FFF2C2",
      "Energy metabolism" = "#B5E6DE",
      "Synaptic vesicle cycle" = "#D3C9EB",
      "Extracellular matrix" = "white",
      "Hormonal and circadian regualtion" = "#E8DCCF",
      "Other" = "white"
    )
  } else if (ptm == "ace") {
    col_value <- c(
      "RNA splicing" = "white",
      "Cytoskeleton & Myelin" = "#FFF2C2",
      "Energy metabolism" = "#B5E6DE",
      "Synaptic vesicle cycle" = "white",
      "Extracellular matrix" = "white",
      "Hormonal and circadian regualtion" = "white",
      "Other" = "white"
    )
  }

  show_num <- 5
  if (ptm == "ubiq") {
    show_num <- 5
  }
  label_list <- plot_df %>%
    filter(
      pval <= genome_wide_p_val
    ) %>%
    group_by(colors) %>%
    slice_min(pval, n = show_num, with_ties = F) %>%
    ungroup() %>%
    filter(colors %in% names(col_value)) %>%
    pull(rsid)

  label_list <- unique(c(label_list, top_hits$rsid))

  plot_df <- plot_df %>%
    mutate(
      label = if_else(rsid %in% label_list, rsid, "")
    )

  dim(plot_df)
  plot_df$label <- str_replace(plot_df$label, "Tau441", "TAU")
  plot_df$label <- str_replace(plot_df$label, "_", "\n")

  tmp_upper <- plot_df[plot_df$beta > 0, ]
  tmp_lower <- plot_df[plot_df$beta <= 0, ]

  my_upper_colors <- c("#656565", "#bfbfbf")

  upper_plot <- ggplot(
    data = plot_df[which(plot_df$beta > 0), ],
    aes(x = rel_pos, y = logged_p)
  ) +
    geom_point(
      aes(color = as.factor(chr)),
      size = 0.8,
      show.legend = F
    ) +
    scale_color_manual(values = rep(my_upper_colors, nrow(axis_df))) +
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
      yintercept = -log10(genome_wide_p_val),
      color = "red",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    labs(x = "", y = bquote("-log"[10] * "(P)")) +
    theme_classic() +
    geom_label_repel(
      data = tmp_upper,
      aes(x = rel_pos, y = logged_p, label = label, fill = colors),
      size = 4,
      segment.size = 0.2,
      point.padding = 0.3,
      ylim = c(maxp / 3, NA),
      box.padding = .5,
      max.overlaps = Inf,
      force = 5,
      force_pull = 0.05,
      min.segment.length = 0,
      segment.curvature = -0.1,
      segment.ncp = 3,
      seed = 2025,
      show.legend = F,
      nudge_x = 0,
      nudge_y = 0,
      direction = "both"
    ) +
    scale_fill_manual(
      values = col_value
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 6),
      # axis.title.y = element_text(margin = margin(r = 3)),
      plot.margin = margin(t = 10, l = 10, r = 10, b = 0)
    )

  print(upper_plot)

  pdf(
    sprintf("%s_%s_miamiplot_acat_20260114.pdf", ptm, ad.trait),
    width = 10,
    height = 4
  )
  plot(upper_plot)
  dev.off()
}
