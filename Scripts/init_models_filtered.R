############## first script to load models (compl>=90%, cont. =<2) #############
library(MicrobiomeGS2)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")


Metadata <- fread("/mnt/nuuk/2021/HRGM/REPR_Genomes_metadata.tsv")

relGenomes <- Metadata[`Completeness (%)`>= 90 & `Contamination (%)` <=2 & !grepl("^d__Archaea", `GTDB Taxonomy`), `HRGM name`]

models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/", IDs = relGenomes)

