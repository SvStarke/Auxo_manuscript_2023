#### Popgen init data

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
#sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")


# get QC for MAGs
completeness_cuttoff <- 85
contamination_cutoff <- 2

popgen_MAGQC <- fread("/mnt/nuuk/2022/MR_popgen_MGX/atlas/completeness.tsv")
popgen_rel_mags <- popgen_MAGQC[Completeness >= completeness_cuttoff & contamination_cutoff <= contamination_cutoff, `Bin Id`]
### load models

models <- fetch_model_collection("/mnt/nuuk/2022/MR_popgen_MGX/models/",
                                 IDs = popgen_rel_mags)

# Read abundancy table for popgen (not yet normalalised)
popgen_mags_abun <- t(read.table("/mnt/nuuk/2022/MR_popgen_MGX/atlas/median_coverage_genomes.tsv"))
popgen_mags_abun1 <- data.frame(popgen_mags_abun)
popgen_mags_abun1$MAGs <- row.names(popgen_mags_abun1)

##filtering rel abundany table for fitlered models
popgen_rel_mags <- data.frame(popgen_rel_mags)
popgen_mags_abun1 <- merge(popgen_mags_abun1, popgen_rel_mags, by.x="MAGs", by.y="popgen_rel_mags")



# Calculate realtive abundancy table
popgen_relabun <- data.table(as.table(popgen_mags_abun))
setnames(popgen_relabun, c("model","sample","prop"))
popgen_rel_mags <- popgen_MAGQC[Completeness >= completeness_cuttoff & contamination_cutoff <= contamination_cutoff, `Bin Id`]
popgen_relabun <- popgen_relabun[model %in% popgen_rel_mags]
popgen_relabun[, prop := prop/sum(prop), by = sample]
popgen_norm <- dcast(popgen_relabun, sample ~ model)



