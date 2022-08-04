#Predict Auxotrophies

model.auxo <- predict_auxotrophies(models)

Auxotrophie <- data.frame(model.auxo)
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
