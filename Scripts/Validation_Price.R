##Prediction of auxotrophies in Price genomes
library(MicrobiomeGS2)
library(dplyr)
library(ggplot2)

models <- fetch_model_collection("/home/svenja/Price_genomes_GCF/models")
model.auxo <- predict_auxotrophies(models)
Auxotrophie <- data.frame(model.auxo)
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
Auxotrophy

Auxotrophy_2 <- melt(Auxotrophy, id.vars = "Genomes",
                     value.name = "Prototrophy", variable.name = "Compound")

AA <- unique(Auxotrophy_2$Compound)
k <- 1
l <- list()

for(i in AA) {
  tmp_Auxos <- Auxotrophy_2[Compound == i]
  Perc_Proto <- ((nrow(tmp_Auxos[tmp_Auxos$Prototrophy == 1]) / 124) * 100)
  tmp_perc_table <- data.frame(Perc_Proto,i)
  l[[k]] <- tmp_perc_table
  k <- k+1
}
Proto_Perc <- rbindlist(l)
Proto_Perc

#visualization

Price_Vali <- ggplot(Proto_Perc, aes(i, Perc_Proto))+
  geom_bar(stat = "identity") +
  xlab("Amino Acid") +
  ylab("Percentage of prototrophies in genomes [%]")+
  theme_bw()
Price_Vali

ggsave("output/plots/Price_Val.pdf", plot = Price_Vali,
       width = 6,height =4)

#validation of loop
(nrow(Auxotrophy_2[Auxotrophy_2$Prototrophy == 1 & Auxotrophy_2$Compound == "Trp"]))
  