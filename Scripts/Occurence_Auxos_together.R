##################   occurence of auxotrophies together ########################

Auxo <- Auxotrophy 
Auxo$Genomes <- NULL
View(Auxo)

AA <- unique(colnames(Auxo))
AA2 <- unique(colnames(Auxo))
new <- list()
k <- 1
for (i1 in AA) {
  print(i1)
  for (i2 in AA2) {
  nauxo <- table(Auxo[[i1]],Auxo[[i2]])
  occu <- data.frame(nauxo)
  occu$A1 <- i1
  occu$A2 <- i2
  occu$perc <- nrow(Auxo[Auxo[[i1]]== 0])
  new[[k]] <- occu
  k <- k + 1
  }
}

occurence <- rbindlist(new)
occurence2 <- occurence[!(occurence$Var1 == 1 | occurence$Var2 == 1), ]
#überlegen statt Frequency Prozentwert berechnen
occurence2$Perc <- ifelse(occurence2$A1 == occurence2$A2, NA, occurence2$Freq/nrow(Auxo) *100)
occurence2$W <- ifelse(occurence2$A1 == occurence2$A2, NA ,occurence2$Freq/occurence2$perc)
View(occurence2)
####################      visualization      ###################################
View(occurence2)
library(ggplot2)
###Frequencies
ggplot(occurence2, aes(A1,A2, fill = Freq)) +
  geom_tile()

###Frequencies in relation to the number of genomes
ggplot(occurence2, aes(A1,A2, fill = Perc)) +
  geom_tile()

### amino acid auxotrophies occuring together in relation to their sole occurence
ggplot(occurence2, aes(A1,A2, fill = W)) +
  geom_tile() 




