##########  association of auxotrophic bacteria with nutritional data ##########

nutr_info <- fread("/mnt/nuuk/2021/HRGM/2020-031_FOC_nutrintake.csv")

#### analysis for amino acids 
nutr_info_AA <- nutr_info[,c(1:25,70)]

####   AA/day
library(dplyr)
nutr_info_AA %>% mutate(across(c(3:25), .fns = ~./GJ))
