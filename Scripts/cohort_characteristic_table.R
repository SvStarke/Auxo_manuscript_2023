#####        cohort characteristic table               #######

###                  this study                      #########

dzhk_info <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/DZHK_finaler_export_v5_mod.csv")



##stool samples with metagenomes
dzhk_info1 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/Metagenome_DZHK_NGS_EMGE_sampleID_delete_EMGE173.tsv", header = F)
setnames(dzhk_info1, c("EMGE","sample"))
dzhk_info2 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/DZHK_finaler_export_v5_mod.csv")
dzhk_info <- merge(dzhk_info1,dzhk_info2)
describe(dzhk_info2$EMGE)



###age
quantile(dzhk_info$`Alter bei Einschluss (Jahre)`, probs = c(0.25, 0.5, 0.75))

quantile(dzhk_info$`Alter bei Einschluss (Jahre)`)

##BMI
quantile(dzhk_info$BMI, probs = c(0.25, 0.5, 0.75))

#### Gender
female_perc <- nrow(dzhk_info[dzhk_info$Geschlecht == "w"]) / nrow(dzhk_info) *100
round((female_perc), digits = 1) 





