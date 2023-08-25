#########      VALIDATION - Experimental vs AGORA2 vs GAPSEQ     #########


##experimental auxotrophy predictions
experiment_auxo <- fread("/home/svenja/workspace/2022/auxotrophies_hrgm/Experimentally_verified_Auxotrophies.csv")
experiment_auxo2 <- experiment_auxo[,-c(1,2,3)]
experiment_auxo2$`GCF-Number` <- as.character(experiment_auxo2$`GCF-Number`)
experiment_auxo2 <- as.data.table(experiment_auxo2)
experiment_auxo2[26,1] <- "GCF_000332015.1"
experiment_auxo2[15,1] <- "GCF_000195955.2"
experiment_auxo2[30,1] <- "GCF_900537995.1"

Auxotrophy_2_exp <- melt(experiment_auxo2, id.vars = "GCF-Number",
                     value.name = "Prototrophy", variable.name = "Compound")

##gapseq auxotrophy predictions
models_gapseq <- fetch_model_collection("/mnt/nuuk/2023/Validation_isolates/models")

auxo_gapseq <- predict_auxotrophies(models_gapseq)
Auxotrophie <- data.frame(auxo_gapseq)
# head(Auxotrophie) 
# summary(Auxotrophie)
# str(Auxotrophie)
# is.data.frame(Auxotrophie)

#column und rows tauschen
Auxotroph <- t(Auxotrophie)
Auxotroph[which(is.na(Auxotroph), arr.ind = T)] <- 1
#is.matrix(Auxotroph)
#data frame erzeugen
Auxotrophy <- data.frame(Auxotroph)
# is.data.frame(Auxotrophy)
# str(Auxotrophy)
# Auxotrophy
Genome <- rownames(Auxotrophy)
Auxotrophy$Genomes <- Genome
# ----
Auxotrophy <- as.data.table(Auxotrophy)

Auxotrophy$Genomes <- gsub('.{8}$', '', Auxotrophy$Genomes)

Auxotrophy_2 <- melt(Auxotrophy, id.vars = "Genomes",
                     value.name = "Prototrophy", variable.name = "Compound")
fwrite(Auxotrophy, file = "Auxotrophy_Validation_gapseq.csv")

###           validation gapseq vs experimental      ###
AA <- unique(Auxotrophy_2$Compound)
li <- list()
k <- 1
for(AAi in AA){
  tmp_1 <- Auxotrophy_2[Compound == AAi]
  tmp_2 <- Auxotrophy_2_exp[Compound == AAi]
  colnames(tmp_2) <- c("GCF-Number", "Comopound", "Prototrophy_exp")
  new <- merge(tmp_1, tmp_2, by.x="Genomes", by.y="GCF-Number")
  li[[k]] <- new
  k <- k+1
  
}
all_Auxo <- rbindlist(li)
all_Auxo <- all_Auxo[,-c(4)]

cont_table <- table(all_Auxo$Prototrophy_exp, all_Auxo$Prototrophy)


Sensitivty <- cont_table[1,1] / (cont_table[1,1] + cont_table[2,1])
Specificity <- cont_table[2,2] / (cont_table[2,2] + cont_table[1,2])
Accuracy <- (cont_table[1,1] + cont_table[2,2]) / (cont_table[1,1] + cont_table[1,2] + cont_table[2,1] + cont_table[2,2])


##AGORA2 auxotrophy predictions
agora_auxo <- fread("/home/svenja/workspace/2022/auxotrophies_hrgm/AGORA2_verified_Auxotrophies.csv")
agora_auxo <- agora_auxo[,-c(21)]

Auxotrophy_2_AGORA <- melt(agora_auxo, id.vars = "V22",
                           value.name = "Prototrophy", variable.name = "Compound")
Auxotrophy_2_AGORA <- as.data.table(Auxotrophy_2_AGORA)
Auxotrophy_2_AGORA$Compound <- str_to_title(Auxotrophy_2_AGORA$Compound) 

###           validation gapseq vs experimental      ###
AA <- unique(Auxotrophy_2$Compound)
lis <- list()
k <- 1
for(AAi in AA){
  tmp_1 <- Auxotrophy_2_AGORA[Compound == AAi]
  colnames(tmp_1) <- c("GCF-Number", "Comopound", "Prototrophy_AGORA")
  tmp_2 <- Auxotrophy_2_exp[Compound == AAi]
  colnames(tmp_2) <- c("GCF-Number", "Comopound", "Prototrophy_exp")
  new <- merge(tmp_1, tmp_2, by.x="GCF-Number", by.y="GCF-Number")
  lis[[k]] <- new
  k <- k+1
  
}
all_Auxo2 <- rbindlist(lis)
all_Auxo2 <- all_Auxo2[,-c(4)]

cont_table2 <- table(all_Auxo2$Prototrophy_exp, all_Auxo2$Prototrophy_AGORA)


Sensitivty <- cont_table2[1,1] / (cont_table2[1,1] + cont_table2[2,1])
Specificity <- cont_table2[2,2] / (cont_table2[2,2] + cont_table2[1,2])
Accuracy <- (cont_table2[1,1] + cont_table2[2,2]) / (cont_table2[1,1] + cont_table2[1,2] + cont_table2[2,1] + cont_table2[2,2])

