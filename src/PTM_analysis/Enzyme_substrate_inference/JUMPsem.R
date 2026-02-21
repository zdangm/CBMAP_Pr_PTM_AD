suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(openxlsx)
  library(magrittr)
  library(patchwork)
  library(limma)
  library(UpSetR)
  library(gridExtra)
  library(qs2)
  library(optparse)
  library(clusterProfiler)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(pheatmap)
  library(impute)
})
rm(list = ls())
conflicted::conflict_prefer_all("dplyr")
NUM_COV <- 56
source("helper_func.R")


enz.ls <- list(
  phos = "phos_JUMPsem.txt",
  ubiq = "ubiq_JUMPsem.txt",
  ace = "ace_JUMPsem.txt"
)


ptm <- "phos"
ptm <- "ubiq"
ptm <- "ace"


enz <- fread(enz.ls[[ptm]]) %>% as.data.frame()
enz <- enz %>% arrange(ADNC_LMH)
rownames(enz) <- enz$id
head(enz[, 1:10])
NUM_COV <- 56
df.ADNC <- enz[, 1:(ncol(enz) - NUM_COV)]


df.ADNC.info <- enz[, (ncol(enz) - NUM_COV + 1):ncol(enz)]


sig_enz <- NULL

phenoData <- enz[, (ncol(enz) - NUM_COV + 1):ncol(enz)]
exprsMatrix <- enz[, 1:(ncol(enz) - NUM_COV)]
# exprsMatrix <- lapply(exprsMatrix, function(ez) {
#   qqnorm(ez, plot.it = F)$x
# })

exprsMatrix <- lapply(exprsMatrix, function(ez) {
  scale(ez)
})

exprsMatrix %<>% as.data.frame() %>% as.matrix() %>% t()

nrow(exprsMatrix)


phe <- "ADNC_LMH"
if (F) {
  phe <- "A_beta_0_3"
  phe <- "Braak_NFT_stage"
  phe <- "B_NFT_braak_0_3"
}
frml <- as.formula(paste0("~ ", phe, " + age + sex_male + PMD + RIN"))
if (F) {
  frml <- as.formula(paste0("~ ", phe))
}

design <- model.matrix(frml, data = phenoData)

colnames(design)
library(limma)
fit <- lmFit(exprsMatrix, design)
fit <- eBayes(fit)
res <- topTable(fit, coef = phe, number = Inf)

res %>%
  rownames_to_column("id") %>%
  write_xlsx(sprintf("../Results/CBMAP_PanNDA/%s_jumpsem_res.xlsx", ptm))

res %>% filter(adj.P.Val < 0.05)
res %>% filter(adj.P.Val < 0.05) %>% nrow()
res %>% filter(P.Value < 0.05)
res %>% filter(P.Value < 0.05) %>% nrow()
ptm
phe
sig_enz <- res %>% filter(adj.P.Val < 0.05) %>% rownames()
fdr_p <- 0.05
if (length(sig_enz) > 0) {
  fdr_p <- res %>% filter(adj.P.Val < 0.05) %>% pull(P.Value) %>% max()
}
maxp <- min(res$P.Value)
DEG <- res
# DEG$regulate<-ifelse(DEG$P.Value>0.05,"unchanged", ifelse(DEG$logFC>0,"up-regulated", ifelse(DEG$logFC<0,"down-regulated","unchanged")))
DEG$regulate <- ifelse(
  DEG$adj.P.Val > 0.05,
  "unchanged",
  ifelse(
    DEG$logFC > 0,
    "up-regulated",
    ifelse(DEG$logFC < 0, "down-regulated", "unchanged")
  )
)
DEG$Genes <- rownames(DEG)


label_df <- DEG %>%
  filter(P.Value <= fdr_p) %>%
  arrange(P.Value) %>%
  slice_head(n = 10)
library(ggrepel)
library(ggthemes)
library(ggbreak)
p <- ggplot(DEG, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(
    aes(color = regulate),
    alpha = 0.5,
    size = 2
  ) +
  geom_hline(
    yintercept = -log10(fdr_p),
    lty = 3,
    col = "black",
    lwd = 0.8
  ) +
  geom_vline(
    xintercept = c(-0.1, 0.1),
    lty = 3,
    col = "black",
    lwd = 0.8
  )
if (ptm != "ubiq") {
  p <- p +
    geom_text_repel(
      data = label_df,
      aes(label = Genes),
      size = 4,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.3,
      force = 19, #
      force_pull = .4, #
      min.segment.length = 0, #
      segment.color = "grey50"
    ) +
    scale_color_manual(
      values = c(
        "down-regulated" = "blue",
        "unchanged" = "grey",
        "up-regulated" = "red"
      )
    ) +
    scale_y_continuous(
      name = bquote("-log"[10] * "(P)"),
      limits = c(0, -log10(maxp) + 1)
    ) +
    theme_few() +
    theme(
      legend.position = "None"
    )
  print(p)
} else {
  p <- p +
    scale_color_manual(
      values = c(
        "down-regulated" = "blue",
        "unchanged" = "grey",
        "up-regulated" = "red"
      )
    ) +
    scale_y_continuous(
      name = bquote("-log"[10] * "(P)"),
      limits = c(0, -log10(maxp) + 1)
    ) +
    theme_few() +
    theme(
      legend.position = "None"
    ) +
    scale_y_break(
      breaks = c(6, 23),
      scales = 0.3,
      ticklabels = c(0, 2, 4, 6, 8, 24, 26, 28),
      space = 0.2
    ) +
    theme(
      axis.ticks.y.right = element_blank(),
      axis.text.y.right = element_blank()
    )
  print(p)
}
print(p)
ggsave(
  filename = sprintf("../Results/CBMAP_PanNDA/JUMPsem_0104_%s.pdf", ptm),
  plot = p,
  width = 4,
  height = 4
)
