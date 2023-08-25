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

# load HRGM model data
Metadata <- fread("data/REPR_Genomes_metadata.tsv")

# abundance data
dzhk_hrgm_abun <- t(read.table("data/mgx_abundances/dzhk_abun.csv"))

rel_models <- intersect(names(which(apply(dzhk_hrgm_abun,1, function(x) any(x > 0)))),
                        Metadata[`Completeness (%)`>= completeness_cuttoff & `Contamination (%)` <= contamination_cutoff, `HRGM name`])

dzhk_hrgm_abun <- dzhk_hrgm_abun[rel_models,]

dzhk_diversity <- data.table(sample = colnames(dzhk_hrgm_abun),
                             D.Shannon = diversity(dzhk_hrgm_abun, MARGIN = 2),
                             D.Simpson = diversity(dzhk_hrgm_abun, MARGIN = 2, index = "simpson"),
                             D.invSimpson = diversity(dzhk_hrgm_abun, MARGIN = 2, index = "invsimpson"),
                             D.richness = specnumber(dzhk_hrgm_abun, MARGIN = 2))


## Calculate realtive abundancy table
dzhk_relabun <- data.table(as.table(dzhk_hrgm_abun))
setnames(dzhk_relabun, c("model","sample","prop"))
dzhk_relabun <- dzhk_relabun[model %in% rel_models]
dzhk_relabun[, prop := prop/sum(prop), by = sample]

# load models
models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/", # /mnt/nuuk/2022/DZHK_MGX/models/
                                 IDs = rel_models)

# Load sample meta information
dzhk_info1 <- fread("data/meta/Metagenome_DZHK_NGS_EMGE_sampleID_delete_EMGE173.tsv", header = F)
setnames(dzhk_info1, c("EMGE","sample"))
dzhk_info2 <- fread("data/meta/DZHK_finaler_export_v5_mod.tsv")
dzhk_info <- merge(dzhk_info1,dzhk_info2)
#describe(dzhk_info2$EMGE)
