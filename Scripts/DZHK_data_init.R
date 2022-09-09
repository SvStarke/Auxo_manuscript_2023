library(vegan)
library(stringr)
library(MicrobiomeGS2)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(MetBrewer)
library(tidyverse)
library(rstatix)
library(RaschSampler)
library(Hmisc)
library(PResiduals)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

# get QC for MAGs
completeness_cuttoff <- 85
contamination_cutoff <- 2

dzhk_MAGQC <- fread("/mnt/nuuk/2022/DZHK_MGX/atlas/completeness.tsv")
dzhk_rel_mags <- dzhk_MAGQC[Completeness >= completeness_cuttoff & contamination_cutoff <= contamination_cutoff, `Bin Id`]
###models

models <- fetch_model_collection("/mnt/nuuk/2022/DZHK_MGX/models/",
                                      IDs = dzhk_rel_mags)

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
dzhk_relabun <- dzhk_relabun[model %in% dzhk_rel_mags]
dzhk_relabun[, prop := prop/sum(prop), by = sample]


# Load sample meta information
dzhk_info1 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/Metagenome_DZHK_NGS_EMGE_sampleID_delete_EMGE173.tsv", header = F)
setnames(dzhk_info1, c("EMGE","sample"))
dzhk_info2 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/DZHK_finaler_export_v5_mod.csv")
dzhk_info <- merge(dzhk_info1,dzhk_info2)


