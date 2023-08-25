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
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")


# get QC for genomes
completeness_cuttoff <- 85
contamination_cutoff <- 2

# load HRGM model data
Metadata <- fread("data/REPR_Genomes_metadata.tsv")

#popgen_MAGQC <- fread("/mnt/nuuk/2022/MR_popgen_MGX/atlas/completeness.tsv")
#popgen_rel_mags <- popgen_MAGQC[Completeness >= completeness_cuttoff & contamination_cutoff <= contamination_cutoff, `Bin Id`]
### load models

# Read abundancy table for popgen (not yet normalalised)
popgen_hrgm_abun <- t(read.table("data/mgx_abundances/troci_abun.csv"))
popgen_hrgm_abun1 <- data.frame(popgen_hrgm_abun)
popgen_hrgm_abun1$MAGs <- row.names(popgen_hrgm_abun1)

rel_models <- intersect(names(which(apply(popgen_hrgm_abun,1, function(x) any(x > 0)))),
                        Metadata[`Completeness (%)`>= completeness_cuttoff & `Contamination (%)` <= contamination_cutoff, `HRGM name`])

popgen_hrgm_abun <- popgen_hrgm_abun[rel_models,]
popgen_hrgm_abun1 <- popgen_hrgm_abun1[popgen_hrgm_abun1$MAGs %in% rel_models,]


models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/",
                                 IDs = rel_models)



##filtering rel abundancy table for fitlered models
rel_models <- data.frame(rel_models)
popgen_hrgm_abun1 <- merge(popgen_hrgm_abun1, rel_models, by.x="MAGs", by.y="rel_models")



# Calculate relative abundancy table
popgen_relabun <- data.table(as.table(popgen_hrgm_abun))
setnames(popgen_relabun, c("model","sample","prop"))
popgen_relabun <- popgen_relabun[model %in% rel_models$rel_models]
popgen_relabun[, prop := prop/sum(prop), by = sample]
popgen_norm <- dcast(popgen_relabun, sample ~ model, value.var = "prop")

