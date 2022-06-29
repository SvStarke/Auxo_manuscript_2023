library(data.table)
library(vegan)
library(stringr)

# Read bundancy table for DZHK (not yet normalalised)
dzhk_mags_abun <- t(read.table("/mnt/nuuk/2022/DZHK_MGX/atlas/median_coverage_genomes.tsv"))

# Calculate alpha-Diversity metrics
dzhk_diversity <- data.table(sample = colnames(dzhk_mags_abun),
                             D.Shannon = diversity(dzhk_mags_abun, MARGIN = 2),
                             D.Simpson = diversity(dzhk_mags_abun, MARGIN = 2, index = "simpson"),
                             D.invSimpson = diversity(dzhk_mags_abun, MARGIN = 2, index = "invsimpson"),
                             D.richness = specnumber(dzhk_mags_abun, MARGIN = 2))

# Calculate realtive abundancy table
dzhk_relabun <- data.table(as.table(dzhk_mags_abun))
setnames(dzhk_relabun, c("model","sample","prop"))
dzhk_relabun[, prop := prop/sum(prop), by = sample]


# Load sample meta information
dzhk_info1 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/Metagenome_DZHK_NGS_EMGE_sampleID_delete_EMGE173.tsv", header = F)
setnames(dzhk_info1, c("EMGE","sample"))
dzhk_info2 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/DZHK_finaler_export_v5_mod.csv")
dzhk_info <- merge(dzhk_info1,dzhk_info2)
rm(dzhk_info1)
rm(dzhk_info2)
