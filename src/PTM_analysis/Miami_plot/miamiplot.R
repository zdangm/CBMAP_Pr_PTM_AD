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
  library(patchwork)
})
rm(list = ls())
conflicted::conflict_prefer_all("dplyr")
ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  mirror = "useast"
)


ptm <- "phospho"
ptm <- "ubiq"
ptm <- "ace"
RES_SUB_DIR <- paste0("../Results/", ptm)
df <- read_xlsx(paste0(
  RES_SUB_DIR,
  "/All.sites.All.traits.intensity.limma.xlsx"
))
colnames(df)[1] <- "identifier"

df$loc_on_MAPT_8 <- if_else(is.na(df$loc_on_MAPT_8), "NA", df$loc_on_MAPT_8)
df$loc_on_MAPT_8 <- if_else(
  str_detect(df$loc_on_MAPT_8, "^not"),
  "NA",
  df$loc_on_MAPT_8
)
df$for_volcano <- if_else(
  df$loc_on_MAPT_8 == "NA",
  df$for_volcano,
  df$loc_on_MAPT_8
)


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

must_cols <- c("identifier", "Gene name", "loc_on_MAPT_8", "for_volcano")
ad.trait <- "lmh"
ad.trait <- "adj.lmh"
ad.trait <- "a3"
ad.trait <- "b3"

df <- df.bak
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
    dplyr::select(all_of(must_cols), starts_with("adj3.ADNC.LMH")) %>%
    rename(
      ADNC.LMH.num.logFC = adj3.ADNC.LMH.num.logFC,
      ADNC.LMH.num.P.Value = adj3.ADNC.LMH.num.P.Value,
      ADNC.LMH.num.adj.P.Val = adj3.ADNC.LMH.num.adj.P.Val
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


{
  genes <- df$`Gene name` %>% unique()
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

  df1 <- df %>% inner_join(gene_info, by = c("Gene name" = "gene_name"))
  dim(df1)
  df1 <- as.data.frame(df1)
  if (F) {
    setdiff(df$`Gene name`, df1$`Gene name`)
    unano_genes <-
      df %>%
      dplyr::filter(identifier %in% setdiff(df$identifier, df1$identifier)) %>%
      arrange(ADNC.LMH.num.adj.P.Val) %>%
      pull(`Gene name`) %>%
      unique()
    df1 %>% filter(`Gene name` %in% unano_genes)
  }

  df1$pos <- str_split_i(df1$for_volcano, "_", 2) %>% as.numeric()
  df1$chr <- factor(df1$chr, levels = c(as.character(1:22), "X", "Y", "MT"))
  df1 <- df1 %>% arrange(chr, start_position, pos)
  df1 <- df1 %>%
    group_by(chr) %>%
    mutate(POS = row_number()) %>%
    ungroup()
  df1 <- df1 %>% dplyr::select(-pos)
  genome_wide_p_val <- df1 %>%
    filter(ADNC.LMH.num.adj.P.Val <= 0.05) %>%
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

  upper_labels <- plot_df %>%
    filter(beta > 0, pval <= genome_wide_p_val) %>%
    arrange(pval)
  message(nrow(upper_labels), " up")
  lower_labels <- plot_df %>%
    filter(beta <= 0, pval <= genome_wide_p_val) %>%
    arrange(pval)
  message(nrow(lower_labels), " down")
  label_list <- rbind(upper_labels, lower_labels)
  nrow(label_list)

  useTOP <- T

  if (useTOP) {
    message("use top hits")

    tmp1 <- plot_df %>%
      filter(beta > 0, pval <= genome_wide_p_val) %>%
      arrange(pval) %>%
      slice_head(n = 15)

    tmp2 <- plot_df %>%
      filter(beta <= 0, pval <= genome_wide_p_val) %>%
      arrange(pval) %>%
      slice_head(n = 15)
    top_hits <- rbind(tmp2, tmp1)
    # label_list <- label_list[!str_detect(label_list$rsid, "^VASP"), ]
  }
  message(sprintf("%d top hits", nrow(top_hits)))

  fine_tuning.cat <- NULL

  if (ptm == "phospho") {
    fine_tuning.cat <- read_excel(
      "../Other_source_data/CBMAP_sigPrs_category_tuning.xlsx",
      sheet = "ph"
    )
  } else if (ptm == "ubiq") {
    fine_tuning.cat <- read_excel(
      "../Other_source_data/CBMAP_sigPrs_category_tuning.xlsx",
      sheet = "ub"
    )
  } else if (ptm == "ace") {
    fine_tuning.cat <- read_excel(
      "../Other_source_data/CBMAP_sigPrs_category_tuning.xlsx",
      sheet = "ac"
    )
  }

  plot_df <- plot_df %>%
    mutate(gene = str_split_i(rsid, "_", 1)) %>%
    left_join(
      fine_tuning.cat,
      by = "gene"
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

  cats_not_white <- col_value[col_value != "white"] %>% names()
  show_num <- 10
  if (ptm == "ubiq") {
    show_num <- 6
  }
  label_list <- plot_df %>%
    filter(
      pval <= genome_wide_p_val
    ) %>%
    group_by(colors) %>%
    slice_min(pval, n = show_num, with_ties = F) %>%
    ungroup() %>%
    filter(colors %in% cats_not_white) %>%
    pull(rsid)

  label_list <- unique(c(label_list, top_hits$rsid))

  plot_df <- plot_df %>%
    mutate(
      label = if_else(rsid %in% label_list, rsid, "")
    )

  plot_df$label <- if_else(
    plot_df$label == "",
    "",
    paste0(
      str_split_i(plot_df$label, "_", 1),
      "_",
      str_split_i(plot_df$label, "_", 3),
      str_split_i(plot_df$label, "_", 2)
    )
  )
  dim(plot_df)
  plot_df$label <- str_replace(plot_df$label, "Tau441", "TAU")
  plot_df$label <- str_replace(plot_df$label, "_", "\n")

  tmp_upper <- plot_df[plot_df$beta > 0, ]
  tmp_lower <- plot_df[plot_df$beta <= 0, ]
  tmp_upper$x <- 0
  tmp_upper$y <- 0
  tmp_lower$x <- 0
  tmp_lower$y <- 0
  upper_nudge <- data.frame(
    x = rep(0, nrow(tmp_upper)),
    y = rep(0, nrow(tmp_upper))
  )
  lower_nudge <- data.frame(
    x = rep(0, nrow(tmp_lower)),
    y = rep(0, nrow(tmp_lower))
  )

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
    labs(x = "", y = bquote("Positive -log"[10] * "(P)")) +
    theme_classic() +
    geom_label_repel(
      data = tmp_upper,
      aes(x = rel_pos, y = logged_p, label = label, fill = colors),
      size = 2.2, #
      segment.size = 0.2,
      point.padding = 0.3,
      ylim = c(maxp / 3, NA),
      box.padding = .5,
      max.overlaps = Inf, #

      force = 19, #
      force_pull = 1.8, #
      min.segment.length = 0, #
      segment.curvature = -0.1, #
      segment.ncp = 3, #
      seed = 2025,
      show.legend = F,
      # nudge_x = -5,
      nudge_y = 2,
      direction = "both" #
    ) +
    scale_fill_manual(
      values = col_value
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      # axis.title.y = element_text(margin = margin(r = 3)),
      plot.margin = margin(t = 10, l = 10, r = 10, b = 0)
    )

  lower_plot <- ggplot(
    data = plot_df[which(plot_df$beta <= 0), ],
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
      position = "top",
      expand = expansion(mult = 0.01)
    ) +
    scale_y_reverse(
      limits = c(maxp * 1.2, 0),
      expand = expansion(mult = c(0.0, 0.02))
    ) +
    geom_hline(
      yintercept = -log10(genome_wide_p_val),
      color = "red",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    labs(x = "", y = bquote("Negative -log"[10] * "(P)")) +
    theme_classic() +
    geom_label_repel(
      data = tmp_lower,
      aes(x = rel_pos, y = logged_p, label = label, fill = colors),
      size = 2.2, #
      segment.size = 0.2,
      point.padding = 0.3,
      ylim = c(NA, -(maxp / 3)),
      box.padding = 0.5,
      max.overlaps = Inf,
      force = 10,
      force_pull = 5,
      min.segment.length = 0.,
      segment.curvature = -0.1,
      segment.ncp = 3,
      seed = 2024,
      direction = "both"
    ) +
    scale_fill_manual(
      values = col_value,
      na.translate = FALSE,
      name = NULL,
      drop = F
    ) +
    theme(
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      # axis.title.y = element_text(margin = margin(r = 3)),
      plot.margin = margin(t = 0, l = 10, r = 10, b = 10)
    )

  #   p <- gridExtra::grid.arrange(
  #     upper_plot,
  #     lower_plot,
  #     nrow = 2
  #   )
  p <- upper_plot / lower_plot + plot_layout(guides = "collect")
  #   p <- p & theme(legend.position = "bottom", legend.box = "horizontal")
  print(p)
  pdf(
    sprintf("%s_%s_miamiplot_batch_20260114.pdf", ptm, ad.trait),
    width = 11,
    height = 5
  )
  plot(p)
  dev.off()
}
