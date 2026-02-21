library(openxlsx)
library(mygene)
library(limma)
library(writexl)
library(biomaRt)

source('/share/home/sunly/Rscript/cbmap_pro/cbmap_586/2_sample_process.R')


# ppi_brain
ppi <- read.xlsx('ppi_list.xlsx')
resa <- queryMany(ppi$protein_A, scopes="uniprot", fields="symbol", species="human")
genea <- data.frame(resa[, c("query", "symbol")])
resb <- queryMany(ppi$protein_B, scopes="uniprot", fields="symbol", species="human")
geneb <- data.frame(resb[, c("query", "symbol")])

ppi$gene_A <- genea$symbol[match(ppi$protein_A, genea$query)]
ppi$gene_B <- geneb$symbol[match(ppi$protein_B, geneb$query)]

ppi_select <- ppi[ppi$gene_A%in%pro_50_raw$gene_name & ppi$gene_B%in%pro_50_raw$gene_name,]
ppi_select <- ppi_select[,c(11,12)]


# ppi_string
ppi_all <- as.data.frame(fread('string_data/9606.protein.physical.links.v12.0.min900.onlyAB.txt'))
ppi_all$protein1 <- sub("9606.", "", ppi_all$protein1)
ppi_all$protein2 <- sub("9606.", "", ppi_all$protein2)
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www")
ensps <- unique(c(ppi_all$protein1, ppi_all$protein2))
id_map <- getBM(
  attributes = c("ensembl_peptide_id", "external_gene_name"),
  filters = "ensembl_peptide_id",
  values = ensps,
  mart = mart
)
colnames(id_map) <- c("protein_id", "gene_symbol")
ppi_all <- merge(ppi_all, id_map, by.x = "protein1", by.y = "protein_id", all.x = TRUE)
colnames(ppi_all)[which(colnames(ppi_all) == "gene_symbol")] <- "gene_A"
ppi_all <- merge(ppi_all, id_map, by.x = "protein2", by.y = "protein_id", all.x = TRUE)
colnames(ppi_all)[which(colnames(ppi_all) == "gene_symbol")] <- "gene_B"


# ppi_combine
ppi_all_select <- ppi_all[ppi_all$gene_A%in%pro_50_raw$gene_name & ppi_all$gene_B%in%pro_50_raw$gene_name,]
ppi_all_select <- ppi_all_select[,c(4,5)]
ppi_all_select <- rbind(ppi_all_select, ppi_select)
ppi_all_select <- unique(ppi_all_select)


# sample select
unique_pro <- union(ppi_all_select$gene_A, ppi_all_select$gene_B)
pro_select_1 <- pro_50_raw[pro_50_raw$gene_name%in%unique_pro,]
pro_select_1 <- pro_select_1[!duplicated(pro_select_1$gene_name),]
rownames(pro_select_1) <- pro_select_1$gene_name

pro_select_1 <- pro_select_1[,colnames(pro_select_1)%in%sample_adnc$id]
pro_select_1_HC <- pro_select_1[,colnames(pro_select_1)%in%sample_adnc$id[sample_adnc$adnc_num==0]]
pro_select_1_LMH <- pro_select_1[,colnames(pro_select_1)%in%sample_adnc$id[sample_adnc$adnc_num%in%c(1,2,3)]]
pro_select_1_MH <- pro_select_1[,colnames(pro_select_1)%in%sample_adnc$id[sample_adnc$adnc_num%in%c(2,3)]]


# cov regression
pro_regress_cov <- pro_select_1
adnc_id <- intersect(sample_adnc$id, colnames(pro_regress_cov))
pro_regress_cov <- pro_regress_cov[,colnames(pro_regress_cov)%in%sample_adnc$id]
pro_regress_cov <- pro_regress_cov[,match(adnc_id, colnames(pro_regress_cov))]

sample_adnc <- sample_adnc[sample_adnc$id %in% adnc_id,]
sample_adnc <- sample_adnc[match(adnc_id, sample_adnc$id),]

design <- model.matrix( ~ adnc_num, data = sample_adnc)
expr_corrected <- removeBatchEffect(pro_regress_cov, 
                                    batch = NULL,
                                    covariates = sample_adnc[, c("age", "sex_male", "PMD", "RIN")],
                                    design = design)

pro_select_2 <- data.frame(expr_corrected)
pro_select_2_HC <- pro_select_2[,colnames(pro_select_2)%in%sample_adnc$id[sample_adnc$adnc_num==0]]
pro_select_2_LMH <- pro_select_2[,colnames(pro_select_2)%in%sample_adnc$id[sample_adnc$adnc_num%in%c(1,2,3)]]
pro_select_2_MH <- pro_select_2[,colnames(pro_select_2)%in%sample_adnc$id[sample_adnc$adnc_num%in%c(2,3)]]


pro_select_1$gene <- rownames(pro_select_1)
pro_select_1 <- pro_select_1[,c(ncol(pro_select_1), 1:(ncol(pro_select_1)-1))]
pro_select_1_HC$gene <- rownames(pro_select_1_HC)
pro_select_1_HC <- pro_select_1_HC[,c(ncol(pro_select_1_HC), 1:(ncol(pro_select_1_HC)-1))]
pro_select_1_LMH$gene <- rownames(pro_select_1_LMH)
pro_select_1_LMH <- pro_select_1_LMH[,c(ncol(pro_select_1_LMH), 1:(ncol(pro_select_1_LMH)-1))]
pro_select_1_MH$gene <- rownames(pro_select_1_MH)
pro_select_1_MH <- pro_select_1_MH[,c(ncol(pro_select_1_MH), 1:(ncol(pro_select_1_MH)-1))]

pro_select_2$gene <- rownames(pro_select_2)
pro_select_2 <- pro_select_2[,c(ncol(pro_select_2), 1:(ncol(pro_select_2)-1))]
pro_select_2_HC$gene <- rownames(pro_select_2_HC)
pro_select_2_HC <- pro_select_2_HC[,c(ncol(pro_select_2_HC), 1:(ncol(pro_select_2_HC)-1))]
pro_select_2_LMH$gene <- rownames(pro_select_2_LMH)
pro_select_2_LMH <- pro_select_2_LMH[,c(ncol(pro_select_2_LMH), 1:(ncol(pro_select_2_LMH)-1))]
pro_select_2_MH$gene <- rownames(pro_select_2_MH)
pro_select_2_MH <- pro_select_2_MH[,c(ncol(pro_select_2_MH), 1:(ncol(pro_select_2_MH)-1))]


write.table(ppi_all_select, 'processed_data/wmds/cbmap_586/physical_9/ppi_select.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_1, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_not_regress_cov_all.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_1_HC, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_not_regress_cov_HC.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_1_LMH, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_not_regress_cov_LMH.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_1_MH, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_not_regress_cov_MH.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_2, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_regress_cov_all.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_2_HC, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_regress_cov_HC.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_2_LMH, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_regress_cov_LMH.txt', sep = '\t', quote = F, row.names = F)
write.table(pro_select_2_MH, 'processed_data/wmds/cbmap_586/physical_9/pro_wmds_regress_cov_MH.txt', sep = '\t', quote = F, row.names = F)


# adnc_limma_lmh
adnc_id <- intersect(sample_adnc$id, colnames(pro_50_raw)[4:ncol(pro_50_raw)])
pro_adnc_limma <- pro_select_1[,colnames(pro_select_1)%in%adnc_id]
pro_adnc_limma <- pro_adnc_limma[, match(adnc_id, colnames(pro_adnc_limma))]
sample_adnc <- sample_adnc[sample_adnc$id%in%adnc_id,]
sample_adnc <- sample_adnc[match(adnc_id, sample_adnc$id), ]


# adnc_limma_lmh(0/1)
design <- model.matrix( ~ adnc + age + as.factor(sex_male) + PMD + RIN, data = sample_adnc)
fit <- lmFit(pro_adnc_limma, design)
fit <- eBayes(fit)
adnc_lim_re <- topTable(fit, coef = "adnc", number = Inf, adjust.method = "BH")
adnc_lim_re$genename <- rownames(adnc_lim_re)
adnc_lim_re$pro <- pro_50_raw$protein[match(adnc_lim_re$genename, pro_50_raw$gene_name)]
nrow(adnc_lim_re[adnc_lim_re$adj.P.Val<0.05,])


# adnc_limma_mh(0/1)
adnc_id_mh <- intersect(sample_adnc$id[sample_adnc$adnc_num%in%c(0,2,3)], colnames(pro_50_raw)[4:ncol(pro_50_raw)])
pro_adnc_limma_mh <- pro_select_1[,colnames(pro_select_1)%in%adnc_id_mh]
pro_adnc_limma_mh <- pro_adnc_limma_mh[, match(adnc_id_mh, colnames(pro_adnc_limma_mh))]
sample_adnc_mh <- sample_adnc[sample_adnc$id%in%adnc_id_mh,]
sample_adnc_mh <- sample_adnc_mh[match(adnc_id_mh, sample_adnc_mh$id), ]

design <- model.matrix( ~ adnc + age + as.factor(sex_male) + PMD + RIN, data = sample_adnc_mh)
fit <- lmFit(pro_adnc_limma_mh, design)
fit <- eBayes(fit)
adnc_lim_re_mh <- topTable(fit, coef = "adnc", number = Inf, adjust.method = "BH")
adnc_lim_re_mh$genename <- rownames(adnc_lim_re_mh)
adnc_lim_re_mh$pro <- pro_50_raw$protein[match(adnc_lim_re_mh$genename, pro_50_raw$gene_name)]
nrow(adnc_lim_re_mh[adnc_lim_re_mh$adj.P.Val<0.05,])

write.table(adnc_lim_re, 'processed_data/wmds/cbmap_586/physical_9/dea_result/physical_9_pro_01_LMH.txt', sep = '\t', quote = F, row.names = F)
write.table(adnc_lim_re_mh, 'processed_data/wmds/cbmap_586/physical_9/dea_result/physical_9_pro_01_MH.txt', sep = '\t', quote = F, row.names = F)


