#Predict Auxotrophies

model.auxo <- lapply(models, FUN = predict_auxotrohies)

Auxotrophie <- data.frame(model.auxo)
head(Auxotrophie) 
summary(Auxotrophie)
str(Auxotrophie)
is.data.frame(Auxotrophie)

#column und rows tauschen
Auxotroph <- t(Auxotrophie)

is.matrix(Auxotroph)
#data frame erzeugen
Auxotrophy <- data.frame(Auxotroph)
is.data.frame(Auxotrophy)
str(Auxotrophy)
Auxotrophy
Genome <- rownames(Auxotrophy)
Auxotrophy$Genomes <- Genome
# ----
Auxotrophy <- as.data.table(Auxotrophy)
