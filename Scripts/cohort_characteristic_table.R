#####        cohort characteristic table               #######

###                  DZHK cohort                      #########

dzhk_info <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/DZHK_finaler_export_v5_mod.csv")

###age
median_Age <- median(dzhk_info$`Alter bei Einschluss (Jahre)`)
round((median_Age), digits = 1)
sd <- sd(dzhk_info$`Alter bei Einschluss (Jahre)`)
round((sd), digits = 1)

##BMI
median_BMI <- median(dzhk_info$BMI)
round((median_BMI), digits = 1)
sd_BMI <- sd(dzhk_info$BMI)
round(sd_BMI, digits = 1)

#### Gender
female_perc <- nrow(dzhk_info[dzhk_info$Geschlecht == "w"]) / nrow(dzhk_info) *100
round((female_perc), digits = 1) 

##stool samples with metagenomes
dzhk_info1 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/Metagenome_DZHK_NGS_EMGE_sampleID_delete_EMGE173.tsv", header = F)
setnames(dzhk_info1, c("EMGE","sample"))
dzhk_info2 <- fread("/mnt/nuuk/2022/DZHK_MGX/sample_info/DZHK_finaler_export_v5_mod.csv")
dzhk_info <- merge(dzhk_info1,dzhk_info2)
describe(dzhk_info2$EMGE)

nrow(dzhk_info)

###                  popgen cohort                  ########


popgen <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/BSPSPC_nondietary_data_rev_date_mod.csv")
##transform to readible format
popgen$sampledate <- format(as.Date(popgen$datum, format = "%d/%m/%Y"), "%Y-%m-%d")
popgen$bday <- format(as.Date(popgen$gebdat, format = "%d/%m/%Y"), "%Y-%m-%d")

### calculate age
library(lubridate)
popgen$age <- time_length(difftime(as.Date(popgen$bday), as.Date(popgen$sampledate)), "years")

popgen2 <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/1.metadata1.tsv")
popgen1 <- merge(popgen2, popgen, by.x="new_id", by.y="new_id")

median_Age <- median(popgen1$age)
round((median_Age), digits = 1)
sd <- sd(popgen1$age)
round((sd), digits = 1)

##BMI
median_BMI <- median(popgen1$BMI)
round((median_BMI), digits = 1)
sd_BMI <- sd(popgen1$BMI)
round(sd_BMI, digits = 1)

#### Gender
popgen2 <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/1.metadata2.tsv")
female_perc <- nrow(popgen2[popgen2$sex == "2"]) / nrow(popgen2)*100
round(print(female_perc), digits = 1) 




