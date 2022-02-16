library(MicrobiomeGS2)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(MetBrewer)
library(tidyverse)
library(rstatix)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

models <- fetch_model_collection("/mnt/nuuk/2021/HRGM/models/")

