library(data.table)

setwd('/data/shared_data/China_Brain_MultiOmics/humanBrain_protein/processed/')


# data input
pro_raw <- as.data.frame(fread('processed_data/human_brain_pro_log_729.txt'))
pro_50_raw <- read.delim("processed_data/human_brain_pro_729_50_norm_rint.txt", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
pro_50_raw$gene_protein <- rownames(pro_50_raw)
pro_50_raw$gene_name <- sub("_.*", "", pro_50_raw$gene_protein)
pro_50_raw$protein <- sub(".*_", "", pro_50_raw$gene_protein)
pro_50_raw <- pro_50_raw[,c(730:732,1:729)]
sample_infor <- read.csv('/data/shared_data/China_Brain_MultiOmics/sample_information/final/CBMAP_sample_info_1187_final_20250523.csv')


# sex mismatch sample
sample_out <- c('A2017CBB026','PTB350','XYA20200730','A2023CBB078','A2019CBB015')
pro_raw <- pro_raw[,-which(colnames(pro_raw)%in%sample_out)]
pro_50_raw <- pro_50_raw[,-which(colnames(pro_50_raw)%in%sample_out)]


# sample information prepare
sample_infor$id <- ifelse(sample_infor$bank=='zju', paste0('A',sample_infor$id), sample_infor$id)
sample_infor <- sample_infor[-which(sample_infor$id%in%sample_out),]
sample_infor$adnc <- ifelse(sample_infor$ADNC==1&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0, 1,
                            ifelse(sample_infor$ADNC==0&sample_infor$LBD==0&sample_infor$ARTAG_Gray.matter==0&sample_infor$ARTAG_Perivascular==0&sample_infor$ARTAG_Subpial==0&sample_infor$OTHER==0&sample_infor$diag_dementia==0&sample_infor$diag_PD==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0, 0, NA))
sample_infor$adnc_num <- ifelse(sample_infor$adnc==1 & sample_infor$ADNC.L==1, 1,
                                ifelse(sample_infor$adnc==1 & sample_infor$ADNC.M==1, 2,
                                       ifelse(sample_infor$adnc==1 & sample_infor$ADNC.H==1, 3, 
                                              ifelse(sample_infor$adnc==0 & sample_infor$ADNC==0, 0, NA))))
sample_infor$part <- ifelse(sample_infor$ADNC==0&sample_infor$PART==1&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0,1,
                            ifelse(sample_infor$ADNC==0&sample_infor$PART==0&sample_infor$ARTAG==0&sample_infor$OTHER==0&sample_infor$diag_dementia==0&sample_infor$diag_PD==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0,0, NA))
sample_infor$late <- ifelse(sample_infor$LATE==1&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0,1,
                            ifelse(sample_infor$LATE==0&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0,0, NA))
sample_infor$artag <- ifelse(sample_infor$ARTAG==1&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0,1,
                             ifelse(sample_infor$PART==0&sample_infor$ARTAG==0&sample_infor$OTHER==0&sample_infor$diag_dementia==0&sample_infor$diag_PD==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0,0, NA))
sample_infor$lbd <- ifelse(sample_infor$LBD==1&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0, 1,
                           ifelse(sample_infor$ADNC==0&sample_infor$LBD==0&sample_infor$ARTAG_Gray.matter==0&sample_infor$ARTAG_Perivascular==0&sample_infor$ARTAG_Subpial==0&sample_infor$OTHER==0&sample_infor$diag_dementia==0&sample_infor$diag_PD==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0, 0, NA))
sample_infor$cvd <- ifelse(sample_infor$CVD==1&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0, 1,
                           ifelse(sample_infor$CVD==0&sample_infor$OTHER==0&sample_infor$diag_ALS==0&sample_infor$diag_SCZ==0&sample_infor$diag_epilepsy==0&sample_infor$diag_others==0, 0, NA))

sample_infor_select <- sample_infor[sample_infor$id%in%colnames(pro_50_raw),]
# sample_infor_select_1 <- sample_infor_select[!is.na(sample_infor_select$ADNC),]


# braak
sample_braak <- sample_infor[!is.na(sample_infor$Braak.NFT.stage),]
sample_braak <- sample_braak[which(sample_braak$OTHER==0&sample_braak$diag_ALS==0&sample_braak$diag_epilepsy==0&sample_braak$diag_others==0&sample_braak$diag_SCZ==0),]
sample_braak$bscore <- ifelse(sample_braak$Braak.NFT.stage==0, 0,
                              ifelse(sample_braak$Braak.NFT.stage >= 1 & sample_braak$Braak.NFT.stage <= 2, 1,
                                     ifelse(sample_braak$Braak.NFT.stage >= 3 & sample_braak$Braak.NFT.stage <= 4, 2, 3)))
sample_braak$braak_fac <- factor(sample_braak$Braak.NFT.stage, ordered = TRUE)


# abeta
sample_abeta <- sample_infor[!is.na(sample_infor$A_beta_0_3),]
sample_abeta <- sample_abeta[which(sample_abeta$OTHER==0&sample_abeta$diag_ALS==0&sample_abeta$diag_epilepsy==0&sample_abeta$diag_others==0&sample_abeta$diag_SCZ==0),]


#cscore
sample_cscore <- sample_infor[sample_infor$bank%in%c('zju','pumc'),]
sample_cscore <- sample_cscore[!is.na(sample_cscore$C_cerad_0_3),]
sample_cscore <- sample_cscore[which(sample_cscore$OTHER==0&sample_cscore$diag_ALS==0&sample_cscore$diag_epilepsy==0&sample_cscore$diag_others==0&sample_cscore$diag_SCZ==0),]


#adnc
sample_adnc <- sample_infor[sample_infor$adnc%in%c(0,1),]
sample_adnc$adnc_num <- ifelse(sample_adnc$adnc==1 & sample_adnc$ADNC.L==1, 1,
                               ifelse(sample_adnc$adnc==1 & sample_adnc$ADNC.M==1, 2,
                                      ifelse(sample_adnc$adnc==1 & sample_adnc$ADNC.H==1, 3, 0)))
sample_adnc$adnc_fac <- factor(sample_adnc$adnc_num, ordered = TRUE) 


# part
sample_part <- sample_infor[sample_infor$part%in%c(0,1),]
sample_part$part <- as.factor(sample_part$part)


# late
sample_late <- sample_infor[sample_infor$late%in%c(0,1),]
sample_late$late <- as.factor(sample_late$late)


# artag
sample_artag <- sample_infor[sample_infor$artag%in%c(0,1),]
sample_artag$artag <- as.factor(sample_artag$artag)


#lbd
sample_lbd <- sample_infor[sample_infor$lbd%in%c(0,1),]


#cvd
sample_cvd <- sample_infor[sample_infor$cvd%in%c(0,1),]





