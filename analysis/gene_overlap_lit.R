# RELEVANT GENES
library(data.table)

setwd("~/Desktop/")

#BCRF

Marina.genes <- read.csv("gene_research/Cand.gene.list.csv", header = F)
GWAS_climate <- read.csv("gene_research/GWAS_climate_genes.csv")
GWAS_beak <- read.csv("gene_research/GWAS_bill_genes.csv")
GWAS_body <- read.csv("gene_research/GWAS_body_genes.csv")
Schweizer_lit <- read.csv("gene_research/Schweizer et al 2021_supplementary_data/Table_S9_simplified.csv")
BCRF_genes <- read.csv("/Users/ericarobertson/Desktop/BCRF/Desaix/BCRF_gene.csv")
alt_genes <- read.csv("/Users/ericarobertson/Downloads/18jan21_mlim_JH_supptables_revised/mlim_JH_TableS2.csv")
Lawson <- read.csv("gene_research/Lawson sup.csv")

Funk.cheek <- read.csv("/Users/ericarobertson/Downloads/cheek.csv")
Funk.crown <- read.csv("/Users/ericarobertson/Downloads/crown.csv")
Funk.body <- read.csv("/Users/ericarobertson/Downloads/body.csv")


fst_genes <- read.csv("GCRF/GCRF_Pub/figures/fst_snp_table_BCRF.csv")
lmm_genes <- read.csv("GCRF/GCRF_Pub/results/bcftools_BCRF/genes_trait_lmm0.5_25kb_BCRF_soft_new_clean.csv")
bslmm_genes <- read.csv("GCRF/GCRF_Pub/results/bcftools_BCRF/genes_trait_bslmm0.5_0.01_25kb_BCRF_soft_new_clean.csv")

overlap_genes <- intersect(BCRF_genes$x, bslmm_genes$gene_name)
overlap_genes <- intersect(Lawson$fortis.scandens.GLM, bslmm_genes$gene_name)
overlap_genes <- intersect(toupper(Schweizer_lit$gene_name_pman), bslmm_genes$gene_name)
overlap_genes <- intersect(Marina.genes$V1, lmm_genes$gene_name)
overlap_genes <- intersect(Funk.cheek$Table.S3.1.Genes.within.10kb.of.SNP.significantly.associated.with.cheek.color., bslmm_genes$gene_name)
overlap_genes <- intersect(Funk.crown$Table.S3.2.Genes.within.10kb.of.SNP.significantly.associated.with.crown.color., bslmm_genes$gene_name)
overlap_genes <- intersect(Funk.body$Table.S3.3.Genes.within.10kb.of.SNP.significantly.associated.with.body., bslmm_genes$gene_name)

cat("Overlapping genes:", overlap_genes, "\n")

# Create an empty list to store results
overlap_results <- list()

# Loop over each column in the Lawson data frame
for (col_name in colnames(Lawson)) {
  overlap_results[[col_name]] <- intersect(Lawson[[col_name]], fst_genes$gene_name)
}

# Print or inspect the results
overlap_results


lmm_genes <- read.csv("GCRF/GCRF_Pub/results/bcftools_BCRF/genes_trait_lmm0.5_25kb_BCRF_soft_new.csv")
bslmm_genes <- read.csv("GCRF/GCRF_Pub/results/bcftools_BCRF/genes_trait_bslmm0.5_0.01_25kb_BCRF_soft_new.csv")
