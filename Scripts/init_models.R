library(MicrobiomeGS2)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggpubr)
library(MetBrewer)
library(tidyverse)
library(rstatix)
library(Hmisc)
library(PResiduals)
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

###HRGM
models <- fetch_model_collection("/mnt/nuuk/2022/HRGM/models/")



