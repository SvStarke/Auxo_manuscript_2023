library(MicrobiomeGS2)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/")