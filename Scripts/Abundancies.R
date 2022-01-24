##############     Abundancies of amino acid auxotrophies     ##################
Auxo_frequency <- table(Auxotrophy_2$Compound, Auxotrophy_2$Prototrophy)
Freq_auxos <- as.data.frame.matrix(Auxo_frequency)
test <- as.data.frame.matrix(Auxo_frequency)
colnames(Freq_auxos) <- c("A", "P")
Freq_auxos$percentage <- (Freq_auxos$A / (nrow(Auxotrophy_2)/21)) * 100


Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_A"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_B"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum== "Firmicutes_G"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_I"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Bdellovibrionota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Cyanobacteria"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Desulfobacterota_A"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Elusimicrobiota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Eremiobacterota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Halobacterota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Myxococcota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Patescibacteria"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Synergistota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Thermoplasmatota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Verrucomicrobiota"] <- "Other"


relAuxos <- unique(Auxotrophy_2$Compound)
relphyla <- c("Actinobacteriota", "Bacteroidota",  "Campylobacterota",
                            "Other", "Firmicutes", "Fusobacteriota", "Proteobacteria")
numbauxos <- list()
k <- 1
for (pi in relphyla) {
  print(pi)
  for(auxoi in relAuxos) {
    x <- Auxotrophy_2[phylum == pi & Compound == auxoi]
    t <- nrow(x[x$Prototrophy == 0])
    n <- data.frame(t)
    n$AA <- auxoi
    n$Phyla <- pi
    numbauxos[[k]] <- n
    k <- k + 1
  }
}

numbauxo <- rbindlist(numbauxos)
numbauxo$perc <- numbauxo$t/ nrow(Auxotrophy)*100

#######################     visualization      #################################
brewer_palette <- c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4")

pt <- ggplot(numbauxo, aes(AA, perc, fill = Phyla)) +
  geom_bar(stat="identity") +
  ylab("Auxotrophies [%]") +
  coord_cartesian(ylim = c(1,100)) +
  xlab("Amino acids") +
  theme(panel.background = element_rect(fill="white", colour= "black")) +
  scale_fill_manual(values = brewer_palette)+
  ggtitle("Abundancies of amino acid auxotrophies in HRGM genomes")
pt
