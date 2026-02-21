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
    library(impute)
    library(ggrepel)
    library(clusterProfiler)
    library(tidyverse)
    library(conflicted)
})
rm(list = ls())
conflict_prefer_all("dplyr")
rep_code_2_id <- data.frame(repeat_code = paste0("repeat", 1:7), id = c("PTB425", "PTB087", "XYA20221024", "XYA20191227", "A2024CBB001", "A2023CBB072", "A2022CBB082"))
rep_code_2_id <- rep_code_2_id %>% filter(!repeat_code %in% c("repeat2", "repeat6"))

setwd("/share/home/lik/Scripts/likun/Random_tasks/Scripts/")
RES_SUB_DIR <- "../Results/CVs_for_reps/"
intensity.ls <- list(
    pr = "/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/XB07045B4DA_0311/XB07045B4DA_mix+sample/MS_identified_information.txt",
    phospho = "/data/shared_data/China_Brain_MultiOmics/humanBrain_Phospho/XB07045B4DPST_0311/XB07045B4DPST_mix+sample/L0G0/MS_identified_information.txt",
    ubiq = "/data/shared_data/China_Brain_MultiOmics/huamnBrain_ubiquitylation/XB07045B4DPUb_初分析报告释放/泛素化修饰组结果/MS_identified_information.txt",
    ace = "/data/shared_data/China_Brain_MultiOmics/huamnBrain_Acetylation/乙酰化修饰结果/MS_identified_information.txt"
)

intensity.ls <- list(
  pr = "/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/XB07045B4DA_0311/XB07045B4DA_mix+sample/MS_identified_information.txt"
)
my_hist_plot <- function(long_df) {
    p <- ggplot(long_df, aes(x = cv)) +
        geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black") +
        theme_bw() 
    p
}
i <- 1
for (i in seq_along(intensity.ls)) {
    ptm <- names(intensity.ls)[i]
    ms_info <- fread(intensity.ls[[i]])
    ms_info <- ms_info %>% select(1:7, starts_with("Intensity"))
    colnames(ms_info) <- gsub("Intensity ", "", colnames(ms_info))
    if (ptm == "pr") {
      ms_info <- ms_info %>% 
          mutate(id = sprintf("%s_%s", `Gene name`, `Protein accession`)) %>% 
          select(id, any_of(rep_code_2_id$repeat_code), any_of(rep_code_2_id$id))
    } else {
      ms_info <- ms_info %>%
          mutate(id = sprintf("%s.%s_%s_%s", `Gene name`, `Protein accession`, `Amino acid`, Position)) %>%
          select(id, any_of(rep_code_2_id$repeat_code), any_of(rep_code_2_id$id))
    }
    ms_cvs <- ms_info %>%
        rowwise() %>%
        summarise(
            id = id,
            cv1 = sd(c(repeat1, PTB425)) / mean(c(repeat1, PTB425)),
            cv3 = sd(c(repeat3, XYA20221024)) / mean(c(repeat3, XYA20221024)),
            cv4 = sd(c(repeat4, XYA20191227)) / mean(c(repeat4, XYA20191227)),
            cv5 = sd(c(repeat5, A2024CBB001)) / mean(c(repeat5, A2024CBB001)),
            cv7 = sd(c(repeat7, A2022CBB082)) / mean(c(repeat7, A2022CBB082))
        )
    log2_ms_info <- cbind(ms_info[, 1], log2(ms_info[, -1] + 1))
    log2_ms_cvs <- log2_ms_info %>%
        rowwise() %>%
        summarise(
            id = id,
            cv1 = sd(c(repeat1, PTB425)) / mean(c(repeat1, PTB425)),
            cv3 = sd(c(repeat3, XYA20221024)) / mean(c(repeat3, XYA20221024)),
            cv4 = sd(c(repeat4, XYA20191227)) / mean(c(repeat4, XYA20191227)),
            cv5 = sd(c(repeat5, A2024CBB001)) / mean(c(repeat5, A2024CBB001)),
            cv7 = sd(c(repeat7, A2022CBB082)) / mean(c(repeat7, A2022CBB082))
        )

    wb <- createWorkbook()

    addWorksheet(wb, "cv")
    addWorksheet(wb, "log2_then_cv")

    writeData(wb, "cv", ms_cvs)
    writeData(wb, "log2_then_cv", log2_ms_cvs)

    saveWorkbook(wb, sprintf("%s%s_cvs_and_logged_cvs.xlsx", RES_SUB_DIR, ptm), overwrite = TRUE)



    p.ms_cvs.ls <- list()
    head(ms_cvs)
    ms_cvs_median <- ms_cvs %>% 
        # filter(!is.na(cv1) & !is.na(cv3) & !is.na(cv4) & !is.na(cv5) & !is.na(cv7)) %>%   # 是否需要去掉那些质量差的位点
        rowwise() %>% 
        summarise(
            id = id, 
            median_cv = median(c_across(starts_with("cv")), na.rm = TRUE)
        )
    head(ms_cvs_median)
    ms_cvs_long <- ms_cvs %>% 
        pivot_longer(-id, names_to = "rep", values_to = "cv") %>% 
        filter(!is.na(cv))

    p.ms_cvs.ls[["all"]] <- my_hist_plot(ms_cvs_long)
    p.ms_cvs.ls[["median"]] <- my_hist_plot(ms_cvs_median %>% select(id, cv = median_cv))
    p.ms_cvs.ls[["cv1"]] <- my_hist_plot(ms_cvs_long %>% filter(rep == "cv1")) 
    p.ms_cvs.ls[["cv3"]] <- my_hist_plot(ms_cvs_long %>% filter(rep == "cv3")) 
    p.ms_cvs.ls[["cv4"]] <- my_hist_plot(ms_cvs_long %>% filter(rep == "cv4")) 
    p.ms_cvs.ls[["cv5"]] <- my_hist_plot(ms_cvs_long %>% filter(rep == "cv5")) 
    p.ms_cvs.ls[["cv7"]] <- my_hist_plot(ms_cvs_long %>% filter(rep == "cv7")) 
    
    marrangeGrob(grobs = p.ms_cvs.ls, ncol = 1, nrow = 1, top = NULL) %>% 
        ggsave(sprintf("%s%s_cvs_hist.pdf", RES_SUB_DIR, ptm), ., width = 6, height = 4)


    log2_ms_cvs_long <- log2_ms_cvs %>% 
        pivot_longer(-id, names_to = "rep", values_to = "cv") %>% 
        filter(!is.na(cv))
    log2_ms_cvs_median <- log2_ms_cvs %>% 
        # filter(!is.na(cv1) & !is.na(cv3) & !is.na(cv4) & !is.na(cv5) & !is.na(cv7)) %>%   # 是否需要去掉那些质量差的位点
        rowwise() %>% 
        summarise(
            id = id, 
            median_cv = median(c_across(starts_with("cv")), na.rm = TRUE)
        )
    p.log2_ms_cvs.ls <- list()
    p.log2_ms_cvs.ls[["all"]] <- my_hist_plot(log2_ms_cvs_long)
    p.log2_ms_cvs.ls[["median"]] <- my_hist_plot(log2_ms_cvs_median %>% select(id, cv = median_cv))
    p.log2_ms_cvs.ls[["cv1"]] <- my_hist_plot(log2_ms_cvs_long %>% filter(rep == "cv1")) 
    p.log2_ms_cvs.ls[["cv3"]] <- my_hist_plot(log2_ms_cvs_long %>% filter(rep == "cv3")) 
    p.log2_ms_cvs.ls[["cv4"]] <- my_hist_plot(log2_ms_cvs_long %>% filter(rep == "cv4")) 
    p.log2_ms_cvs.ls[["cv5"]] <- my_hist_plot(log2_ms_cvs_long %>% filter(rep == "cv5")) 
    p.log2_ms_cvs.ls[["cv7"]] <- my_hist_plot(log2_ms_cvs_long %>% filter(rep == "cv7")) 

    marrangeGrob(grobs = p.log2_ms_cvs.ls, ncol = 1, nrow = 1, top = NULL) %>% 
        ggsave(sprintf("%s%s_log2_then_cvs_hist.pdf", RES_SUB_DIR, ptm), ., width = 6, height = 4)
}

