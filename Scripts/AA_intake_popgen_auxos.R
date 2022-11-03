####### popgen nutritional analysis and frequency of auxotrophic bacteria ######

### read data table of nutritional analysis   ###

nutr_data_popgen <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/11.food.groups.nutrients.f1.f2.tsv")

####calculate energy-%
nutr_data_popgen[,c(6:28)] <- nutr_data_popgen[,c(6:28)] * 17
nutr_data_popgen <- nutr_data_popgen%>% mutate(across(c(6:28), .fns = ~./GJ)*100)


#calculate total amount of every consumed amino acid
ID <- unique(nutr_data_popgen$new_id)
k <- 1
nutr_list <- list()

for(id in ID) {
  tmp <- nutr_data_popgen[new_id == id]
  tmp <- tmp[,c(1,6:28)]
  tmp[is.na(tmp)] <- 0
  tmp_sums <- colSums(tmp[,-1])
  tmp_sums <- data.frame(tmp_sums)
  tmp_sums_t <- data.frame(t(tmp_sums))
  tmp_sums_t$ID <- id
  nutr_list[[k]] <- tmp_sums_t
  k <- k+1
}
nutr_popgen_sums_AA <- rbindlist(nutr_list)



#exclude columns without AA data
nutr_popgen_sums_AA <- nutr_popgen_sums_AA[,!c(5,8,14,15,16)]



sub <- unique(popgen_relabun$sample)
relAA <- unique(Auxotrophy_2$Compound)
Auxotrophy_2

p <- list()
k <- 1


for (subi in sub) {
  print(subi) 
  for (AAi in relAA) {
    x <- popgen_relabun[sample == subi]
    y <- Auxotrophy_2[Compound == AAi]
    z <- merge(x,y, by.x = "model", by.y = "Genomes")
    t <- z[Prototrophy == 0]
    p[[k]] <- t
    k <- k +1
  }
}

u <- rbindlist(p) 

sumfreq <- aggregate(u$prop, by=list(sample=u$sample, AA=u$Compound), FUN=sum)


###prepare data for merging
popgen_F1 <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/1.metadata1.tsv")
popgen_F1$newID <- sub("$", "_F1", popgen_F1$new_id)
#head(popgen_F1)
#colnames(popgen_F1)

popgen_F1 <- popgen_F1[,c(1,2,20)]

##load second timepoint
popgen_F2 <- fread("/mnt/nuuk/2022/MR_popgen_MGX/meta/1.metadata2.tsv")
popgen_F2$newID <- sub("$", "_F2", popgen_F2$new_id)
colnames(popgen_F2) [2] <- "new_id"

popgen_F2 <- popgen_F2[,c(1,2,19)]

F1_F2 <- rbind(popgen_F2, popgen_F1)
#View(popgen_F1)
##merge information about time points with sample info
popgen_samples <- fread("/mnt/nuuk/2022/MR_popgen_MGX/atlas/atlas_samples.csv")
popgen_samples$new_name <- sub("-L001|-L002", "", popgen_samples$Full_Name)
popgen_data <- merge(F1_F2, popgen_samples, by.x= "newID", by.y="new_name")
#View(popgen_data)
#delete sample ID with only  available information about one timepoint in original dataframe about samples
popgen_data <- popgen_data[popgen_data$atlas_name != "S26" & popgen_data$atlas_name != "S38" & popgen_data$atlas_name!= "S66"]
popgen_data

###merge all files together
tmp_samples <- merge(popgen_data, sumfreq, by.x= "atlas_name", by.y="sample")
nutr_popgen_sums_AA <- melt(nutr_popgen_sums_AA, id.vars = "ID",
     value.name = "intake/day", variable.name = "Compound")



##loop for merging of last two files
Comp <- unique(nutr_popgen_sums_AA$Compound)
AA <- unique(tmp_samples$AA)
ID <- unique(nutr_popgen_sums_AA$ID)
List <- list()
k <- 1
for (C in Comp) {
  for(A in AA) {
    for(id in ID) {
    tmp_C <- nutr_popgen_sums_AA[Compound == C]
    tmp_I <- tmp_C[ID == id]
    tmp_A <- tmp_samples[AA == A]
    tmp_CA <- merge(tmp_I, tmp_A, by.x="ID", by.y="newID")
    List[[k]] <- tmp_CA
    k <- k+1
     }
  }
}

Merge_popgen <- rbindlist(List)
View(Merge_popgen)


Merge_popgen[, Time := gsub("^.{9}", "", SampleID)]
Merge_popgen_F1 <- Merge_popgen[Time == "F1"]
Merge_popgen_F2 <- Merge_popgen[Time == "F2"]

#### Spearman correlation analysis F1
Auxo <- unique(Merge_popgen_F1$AA)
Intake <- unique(Merge_popgen_F1$Compound)
spear_list <- list()
k <- 1

for(Aux in Auxo) {
  for(In in Intake) {
    tmp_1 <- Merge_popgen_F1[AA == Aux & Compound == In]
    spear <- cor.test(tmp_1$`intake/day`, tmp_1$x, method = "spearman", exact = FALSE)
    spear_intake <- data.table(Intake_AA = In,
                               Auxo = Aux,
                               rho = spear$estimate,
                               pvalue = spear$p.value)
    spear_list[[k]] <- spear_intake
    k <- k+1
  }
}

spear_Popgen_F1 <- rbindlist(spear_list)


spear_Popgen_F1[, padj := p.adjust(pvalue, method = "fdr")]
spear_Popgen_F1[padj < 0.05, sign.label1 := "P < 0.05"]
##small changes for visualization
spear_Popgen_F1$Intake_AA <- sub('.', '', spear_Popgen_F1$Intake_AA)
spear_Popgen_F1$Intake_AA <- tolower(spear_Popgen_F1$Intake_AA)
spear_Popgen_F1$Intake_AA  <- stringr::str_to_title(spear_Popgen_F1$Intake_AA)


##visualization
nutrition_popgen_F1 <- ggplot(spear_Popgen_F1[Auxo!="Gly"], aes(Auxo, Intake_AA, fill = rho)) +
  geom_tile() +
  geom_point(aes(shape = sign.label1)) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  xlab("Auxotrophy") +
  ylab("Intake of amino acids [E%]") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.text.x = element_text(colour = "black"),
  axis.text.y = element_text(colour = "black"))

ggsave("output/plots/AA_intake_Auxos_popgen_F1.pdf", plot = nutrition_popgen_F1,
       width = 7, height = 5)


#### Spearman correlation analysis F2
Auxo <- unique(Merge_popgen_F2$AA)
Intake <- unique(Merge_popgen_F2$Compound)
spear_list <- list()
k <- 1

for(Aux in Auxo) {
  for(In in Intake) {
    tmp_1 <- Merge_popgen_F2[AA == Aux & Compound == In]
    spear <- cor.test(tmp_1$`intake/day`, tmp_1$x, method = "spearman", exact = FALSE)
    spear_intake <- data.table(Intake_AA = In,
                               Auxo = Aux,
                               rho = spear$estimate,
                               pvalue = spear$p.value)
    spear_list[[k]] <- spear_intake
    k <- k+1
  }
}

spear_Popgen_F2 <- rbindlist(spear_list)


spear_Popgen_F2[, padj := p.adjust(pvalue, method = "fdr")]
spear_Popgen_F2[padj < 0.05, sign.label1 := "P < 0.05"]

##small changes for visualization
spear_Popgen_F2$Intake_AA <- sub('.', '', spear_Popgen_F2$Intake_AA)
spear_Popgen_F2$Intake_AA <- tolower(spear_Popgen_F2$Intake_AA)
spear_Popgen_F2$Intake_AA  <- stringr::str_to_title(spear_Popgen_F2$Intake_AA)

## visualization
nutrition_popgen_F2 <- ggplot(spear_Popgen_F2[Auxo != "Gly"], aes(Auxo, Intake_AA, fill = rho)) +
  geom_tile() +
  geom_point(aes(shape = sign.label1)) +
  scale_fill_gradient2(high = "#ca0020", mid = "white", low = "#0571b0") +
  scale_shape_manual(values = 8, na.translate = FALSE) +
  xlab("Auxotrophy") +
  ylab("Intake of amino acids [E%]") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))


ggsave("output/plots/AA_intake_Auxos_popgen_F2.pdf", plot = nutrition_popgen_F2,
       width = 7, height = 5)


