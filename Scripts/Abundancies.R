##############     Abundancies of amino acid auxotrophies in HRGM ##############

Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_A"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_B"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum== "Firmicutes_G"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Firmicutes_I"] <- "Firmicutes"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Bdellovibrionota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Campylobacterota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Cyanobacteria"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Desulfobacterota_A"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Elusimicrobiota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Eremiobacterota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Fibrobacterota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Halobacterota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Myxococcota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Patescibacteria"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Synergistota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Spirochaeotota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Thermoplasmatota"] <- "Other"
Auxotrophy_2$phylum[Auxotrophy_2$phylum == "Verrucomicrobiota"] <- "Other"

relAuxos <- unique(Auxotrophy_2$Compound)
relphyla <- c("Actinobacteriota", "Bacteroidota",
              "Other", "Firmicutes", "Fusobacteriota", "Proteobacteria")

test <- unique(Auxotrophy_2$phylum)
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
brewer_palette <- c('#feedde','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801','#8c2d04')
library(ggplot2)
pt <- ggplot(numbauxo, aes(AA, perc, fill = Phyla)) +
  geom_bar(stat="identity") +
  ylab("Auxotrophies [%]") +
  xlab("Amino acids") +
  #theme_minimal() +
  theme(legend.position = "bottom") +
  theme(axis.line = element_line(size=0.2)) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.text.x = element_text(size = 10, angle = 45, colour = "black", hjust = 1, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 10, colour = "black"))+
  theme(axis.title.x = element_text(colour = "Black", face = "bold", size = 10, margin = margin(2,0,0,0))) +
  theme(axis.title.y = element_text(colour = "Black", face = "bold", size = 10, margin = margin(0,10,0,0))) +
  theme(legend.title = element_text(colour = "black", size = 10, face = "bold",  vjust = 1, margin = margin(10,5,5,0))) +
  theme(legend.text = element_text(size = 10, colour = "black")) +
  coord_cartesian()+
  scale_fill_manual(values = brewer_palette)
pt

pt1 <- annotate_figure(pt, fig.lab = "A")

ggsave("output/plots/Abundancies_HRGM.pdf", plot = pt1,
       width = 7, height = 5)


###################        optional additional code       ######################
### get abundancies without phylum
Auxo_abundance <- table(Auxotrophy_2$Compound, Auxotrophy_2$Prototrophy)
Freq_auxos <- as.data.frame.matrix(Auxo_abundance)
colnames(Freq_auxos) <- c("A", "P")
Freq_auxos$percentage <- (Freq_auxos$A / (nrow(Auxotrophy_2)/21)) * 100
Freq_auxos$AA <- rownames(Freq_auxos)
Freq_auxos <- data.table(Freq_auxos)
Freq_auxos[AA %in% c("Val","Met","Leu","Ile","Trp","Phe","Lys","His","Thr"), is.essential := "essential"]
Freq_auxos[is.na(is.essential), is.essential := "not essential"]
Freq_auxos[order(Freq_auxos$is.essential)]
Freq_auxos$AA <- factor(Freq_auxos$AA, levels = c("His","Ile","Leu","Lys","Met","Phe","Thr","Trp","Val","Ala",
                                                                      "Arg","Asn","Asp","Chor", "Cys","Gln","Glu","Gly","Pro","Ser","Tyr"))
##visualization
Abun_HRGM <- ggplot(Freq_auxos, aes(AA,percentage, fill = is.essential)) +
   geom_bar(stat = "identity") +
  ylab("Auxotrophies [%]") +
  xlab("Amino acids") +
  guides(fill =guide_legend(title="Essentiality")) +
  #theme_minimal() +
  theme(legend.position = "right") +
  theme(axis.line = element_line(size=0.2)) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.text.x = element_text(size = 16, angle = 45, colour = "black", hjust = 1, margin = margin(10,0,0,0))) +
  theme(axis.text.y = element_text(size = 16, colour = "black"))+
  theme(axis.title.x = element_text(colour = "Black",size = 16, margin = margin(2,0,0,0))) +
  theme(axis.title.y = element_text(colour = "Black",size = 16, margin = margin(0,10,0,0))) +
  theme(legend.title = element_text(colour = "black", size = 16,  vjust = 1, margin = margin(10,5,5,0))) +
  theme(legend.text = element_text(size = 16, colour = "black")) +
  coord_cartesian() +
  scale_fill_manual(values = c("#8c2d04",'#fdae6b'))
Abun_HRGM

ggsave("output/plots/Abundancies_HRGM_18.07.22.pdf", plot = Abun_HRGM,
       width = 10, height = 5)
