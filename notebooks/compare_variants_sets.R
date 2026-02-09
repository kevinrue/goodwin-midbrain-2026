library(VariantAnnotation)
bcftools_variants <- readVcf("workflow-dna-seq-gatk-variant-calling/results/merge/filtered.snps.vcf.gz")
gatk_variants <- readVcf("workflow-dna-seq-gatk-variant-calling/results/genotype_gvcfs/filtered.snps.vcf.gz")
gatk_2025_variants <- readVcf("/ceph/project/cncb/albrecht/snakemake-cncb-genetic-demultiplexing/results/bcftools_merge_filtered_variants/all.vcf.gz")

ranges(bcftools_variants)
ranges(gatk_variants)

co_2025_2026 <- countOverlaps(gatk_2025_variants, gatk_variants)
co_2026_2025 <- countOverlaps(gatk_variants, gatk_2025_variants)
table(co_2025_2026)
table(co_2026_2025)

co_bcf <- countOverlaps(bcftools_variants, gatk_variants)
co_gatk <- countOverlaps(gatk_variants, bcftools_variants)
table(co_bcf)
table(co_gatk)

called(bcftools_variants)

m <- match(bcftools_variants, gatk_variants)
head(m)
table(is.na(m))
