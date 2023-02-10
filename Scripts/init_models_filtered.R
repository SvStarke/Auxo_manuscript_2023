############## first script to load models (compl>=85%, cont. =<2) #############
library(MicrobiomeGS2)
library(stringr)
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


Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")
relGenomes <- Metadata[`Completeness (%)`>= 85 & `Contamination (%)` <=2 & !grepl("^d__Archaea", `GTDB Taxonomy`), `HRGM name`]

models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/", IDs = relGenomes)
