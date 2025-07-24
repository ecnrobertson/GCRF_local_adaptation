#MAKING THE FINAL TABLES

library(tidyverse)
library(dplyr)
library(data.table)

setwd("~/Desktop/GCRF/GCRF_Pub/")

##### GWAS SNPS #####
lmm_snps_fst <- read.csv("results/lmm_results/lmm_snps_0.8_fst_soft.csv")
bslmm_snps_fst <- read.csv("results/bslmm_results/bslmm_snps_0.8_fst_soft.csv")

# lmm_trait_genes <- read.csv("results/bcftools_ZEFI/genes_trait_lmm0.8_25kb_ZEFI_soft.csv")
# bslmm_trait_genes <- read.csv("results/bcftools_ZEFI/genes_trait_bslmm0.8_0.05_25kb_ZEFI_soft.csv")

lmm_trait_genes <- fread("results/bcftools_BCRF/genes_trait_lmm0.8_25kb_BCRF_soft.csv")
bslmm_trait_genes <- fread("results/bcftools_BCRF/genes_trait_bslmm0.8_0.05_25kb_BCRF_soft.csv")

all_bslmm <- read.csv("results/bslmm_results/summary_pip0.05_0.8_soft.csv")
all_lmm <- read.csv("results/lmm_results/lmm_summary_outliers_0.8_soft.csv")

lmm_reorder <- fread("figures/reordered_files/beak_depth0.8_soft.lmm_ZFchr_assoc.txt")
bslmm_reorder <- fread("figures/reordered_files/beak_width0.8.0.05_soft.bslmm_ZFchr_params.txt")

all_lmm <- left_join(all_lmm, lmm_reorder %>% select(ZFCHROM, rs), by="rs")
all_bslmm <- left_join(all_bslmm, bslmm_reorder %>% select(ZFCHROM, rs), by="rs")

lmm_table <- left_join(all_lmm, lmm_trait_genes %>% select(rs, gene_name), by="rs")
lmm_table <- left_join(lmm_table, lmm_snps_fst %>% select(rs, FST), by="rs")
lmm_table <- lmm_table %>% select(trait, rs, ZFCHROM, p_wald, FST, gene_name)
# write.csv(lmm_table, "figures/lmm_snp_table_ZEFI.csv", col.names = F)
write.csv(lmm_table, "figures/lmm_snp_table_BCRF_soft.csv", row.names = F)

bslmm_table <- left_join(all_bslmm, bslmm_trait_genes %>% select(rs, gene_name), by="rs")
bslmm_table <- left_join(bslmm_table, bslmm_snps_fst %>% select(rs, FST), by="rs")
bslmm_table <- bslmm_table %>% select(trait, rs, ZFCHROM, gamma, FST, gene_name)
# write.csv(bslmm_table, "figures/bslmm_snp_table_ZEFI.csv", col.names = F)
write.csv(bslmm_table, "figures/bslmm_snp_table_BCRF_soft.csv", row.names = F)

##### FST SNPS #####
fst_snps <- read.csv("results/fst_OutFLANK/OutFLANK_outliers.22.csv")

fst_reorder <- fread("figures/reordered_files/OutFLANK_fst_results_0.005_ZFchr.txt")

# fst_genes <- read.csv("results/bcftools_ZEFI/genes_fst0.8_25kb_ZEFI.csv")
fst_genes <- read.csv("results/bcftools_BCRF/genes_fst0.8_25kb_BCRF.csv")

all_fst <- left_join(fst_snps, fst_reorder %>% select(ZFCHROM, rs), by="rs")

fst_table <- left_join(all_fst, fst_genes %>% select(rs, gene_name), by="rs")

fst_table <- fst_table %>% select(rs, ZFCHROM, FST, gene_name)
# write.csv(fst_table, "figures/fst_snp_table_ZEFI.csv", col.names = F)
write.csv(fst_table, "figures/fst_snp_table_BCRF.csv", col.names = F)

#####
library(ggplot2)

ggplot(data=bslmm)+
  geom_boxplot(aes(x=as.factor(trait), y=FST))+
  geom_hline(yintercept = 0.00449)

ggplot(data=lmm)+
  geom_boxplot(aes(x=as.factor(trait), y=FST))+
  geom_hline(yintercept = 0.00449)

#####
fst <- read.csv("figures/fst_snp_table_BCRF.csv")
lmm <- read.csv("figures/lmm_snp_table_BCRF_soft.csv") %>% filter(!trait=="tail")
bslmm <- read.csv("figures/bslmm_snp_table_BCRF_soft.csv") %>% filter(!trait=="tail")



###
# GWAS RESULTS
lmm_summary_outliers_0_8_soft <- read_csv("results/lmm_results/lmm_summary_outliers_0.8_soft.csv")
colnames(lmm_summary_outliers_0_8_soft)

lmm_genes <- fread("results/bcftools_BCRF/genes_lmm0.8_25kb_BCRF_soft.txt")

lmm_snps_genes <- left_join(lmm_summary_outliers_0_8_soft, lmm_genes, by="rs") %>% 
  select(chr, rs, ps, p_wald, trait, gene_name)

lmm_snps_fst <- read.csv("results/lmm_results/lmm_snps_0.8_fst_soft.csv")

lmm_table <- left_join(lmm_snps_genes, lmm_snps_fst %>% select(rs, FST), by="rs")

lmm_reorder <- fread("figures/reordered_files/beak_depth0.8_soft.lmm_ZFchr_assoc.txt")

lmm_table_ZF <- left_join(lmm_table, lmm_reorder %>% select(rs, ZFCHROM), by="rs")

lmm_table_final <- lmm_table_ZF %>% select(trait, rs, ZFCHROM, p_wald, FST, gene_name)

write.csv(lmm_table_final, "figures/lmm_snp_table_BCRF_soft_final.csv", row.names = F)


bslmm_summary_pip0.05_0.8_soft <- read_csv("results/bslmm_results/summary_pip0.05_0.8_soft.csv")
colnames(lmm_summary_pip0.05_0.8_soft)

bslmm_genes <- fread("results/bcftools_BCRF/genes_bslmm0.8_25kb_BCRF_soft.txt")

bslmm_snps_genes <- left_join(bslmm_summary_pip0.05_0.8_soft, bslmm_genes, by="rs") %>% 
  select(chr, rs, ps, gamma, trait, gene_name)

bslmm_snps_fst <- read.csv("results/bslmm_results/bslmm_snps_0.8_fst_soft.csv")

bslmm_table <- left_join(bslmm_snps_genes, bslmm_snps_fst %>% select(rs, FST), by="rs")

bslmm_reorder <- fread("figures/reordered_files/beak_depth0.8.0.05_soft.bslmm_ZFchr_params.txt")

bslmm_table_ZF <- left_join(bslmm_table, bslmm_reorder %>% select(rs, ZFCHROM), by="rs")

bslmm_table_final <- bslmm_table_ZF %>% select(trait, rs, ZFCHROM, gamma, FST, gene_name)

write.csv(bslmm_table_final, "figures/bslmm_snp_table_BCRF_soft_final.csv", row.names = F)

# FST RESULTS
fst_snps <- read.csv("results/fst_OutFLANK/OutFLANK_outliers.23FST.csv")

fst_reorder <- fread("figures/reordered_files/OutFLANK_fst_results_0.005_ZFchr.txt")

# fst_genes <- read.csv("results/bcftools_ZEFI/genes_fst0.8_25kb_ZEFI.csv")
fst_genes <- read.csv("results/bcftools_BCRF/genes_0.23fst0.8_25kb_BCRF.csv")

all_fst <- left_join(fst_snps, fst_reorder %>% select(ZFCHROM, rs), by="rs")

fst_table <- left_join(all_fst, fst_genes %>% select(rs, gene_name), by="rs")

fst_table <- fst_table %>% select(rs, ZFCHROM, FST, gene_name)
# write.csv(fst_table, "figures/fst_snp_table_ZEFI.csv", col.names = F)
write.csv(fst_table, "figures/fst_snp_table_BCRF.csv", col.names = F)

#FINAL FINAL TABLE
colnames(lmm_table_final)
colnames(bslmm_table_final)
colnames(fst_table)

lmm_table_final <- read.csv("figures/lmm_snp_table_BCRF_soft_final.csv")
bslmm_table_final <- read.csv("figures/bslmm_snp_table_BCRF_soft_final.csv")
fst_table <- read.csv("figures/fst_snp_table_BCRF.csv")   

lmm_genes <- unique(not_tail$gene_name)
write_lines(lmm_genes, "lmm_genes.txt")

bslmm_genes <- unique(bslmm_table_final$gene_name)
write_lines(bslmm_genes, "bslmm_genes.txt")

fst_genes <- unique(fst_table$gene_name)
write_lines(fst_genes, "fst_genes.txt")
