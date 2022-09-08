#### BMI effect on auxotrophic gut bacteria
auxos_popgen <- merge(Auxotrophy, popgen_relabun, by.x= "Genomes", by.y="model")
popgen_all <- merge(popgen_data, auxos_popgen, by.x="atlas_name", by.y="sample")

sub <- unique(popgen_all$atlas_name)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2
h <- list()
k <- 1

for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- popgen_all[atlas_name == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "Genomes", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    h[[k]] <- t
    k <- k +1
  }
}

info_auxo_popgen <- rbindlist(h) 
info_auxo_popgen <- info_auxo_popgen[Compound != "Gly"]


###load information about F1
popgen_F1 <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/1.metadata1.tsv")
popgen_F1$newID <- sub("$", "_F1", popgen_F1$new_id)
F1_names <- unique(popgen_F1$SampleID)
info_auxo_popgen_F1 <- info_auxo_popgen[info_auxo_popgen$SampleID %in% F1_names]
info_auxo_popgen_F1 <- merge(info_auxo_popgen_F1, popgen_F1, by.x="new_id", by.y="new_id")

###sum up per sample
sumfreq_popgen <- aggregate(info_auxo_popgen_F1$prop, by = list(subject = info_auxo_popgen_F1$atlas_name, AA=info_auxo_popgen_F1$Compound,
                                                           BMI=info_auxo_popgen_F1$BMI,sex=info_auxo_DZHK$Geschlecht,age=info_auxo_DZHK$Age), FUN = sum)
###impact BMI on auxotrophic bacteria
###distribution

hist_BMI_popgen <- hist(popgen_F1$BMI)

ggsave()
sumfreq_popgen <- data.table(sumfreq_popgen)
relAA1 <- unique(sumfreq_popgen$AA)
l <- list()
k <- 1
for (AAi in relAA1){
  print(AAi)
  sumfreq_popgen_BMI <- sumfreq_popgen[AA == AAi]
  spearman_popgen_BMI <- cor.test(sumfreq_popgen_BMI$BMI, sumfreq_popgen_BMI$x, method = "spearman", exact =FALSE)
  t <- data.table(Pvalue = spearman_popgen_BMI$p.value,
                  Corr = spearman_popgen_BMI$estimate,
                  AA = AAi)
  l[[k]] <- t
  k <- k +1
  
}

BMI_popgen_Spear <- rbindlist(l)
BMI_popgen_Spear$padjust = p.adjust(BMI_popgen_Spear$Pvalue, method = "fdr")
BMI_popgen_Spear[padjust < 0.05, sign.label1 := "P < 0.05"]

#
ggplot(sumfreq_popgen, aes(x=AA, y=x)) +
  geom_boxplot()
