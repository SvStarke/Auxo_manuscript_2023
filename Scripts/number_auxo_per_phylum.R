######################      numbers of auxotrophies per phylum #################

count <- unique(Auxotrophy_12$count)
phylum <- unique(Auxotrophy_12$phylum)
Auxotrophy_13 <- Auxotrophy_12
tables <- list()
k <- 1
for (i in count) {
  print(i) 
    for (pi in phylum) {
      x <- Auxotrophy_12[count == i & phylum == pi]
      table <- nrow(x)
      tablea <- data.frame(table)
      tablea$phylum <- pi
      tablea$count <- i
      tablea$Perc <- nrow(Auxotrophy_13[Auxotrophy_13$phylum == pi])
      tables[[k]] <- tablea
      k <- k + 1
  }
}
numb <- rbindlist(tables)
numb$abun <- numb$table / numb$Perc * 100

View(numb)

'#feedde','#fdd0a2','#fdae6b','#fd8d3c','#f16913','#d94801'
#Firmicutes
fi <- ggplot(numb[phylum == "Firmicutes"], aes (count,abun, fill = "#fdae6b")) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(face = "bold")) +
  scale_fill_manual(values = "#fdae6b") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
fi

#Actinobacteria
a <- ggplot(numb[phylum == "Actinobacteriota"], aes (count, abun, fill = "#feedde")) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") + 
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) +
  scale_fill_manual(values = "#feedde") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
a

#Bacteroidetes
b <- ggplot(numb[phylum == "Bacteroidota"], aes (count,abun, fill = "#fdd0a2")) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#fdd0a2") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
b

#Proteobacteria
p <- ggplot(numb[phylum == "Proteobacteria"], aes (count, abun, fill = "#f16913")) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) +
  scale_fill_manual(values = "#f16913") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
p
#Fusobacteriota
fu <- ggplot(numb[phylum == "Fusobacteriota"], aes (count, abun, fill = "#fd8d3c")) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#fd8d3c") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
fu
#Other
ot <- ggplot(numb[phylum == "Other"], aes (count, abun, fill = "#d94801")) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(axis.title.x = element_text(face = "bold")) +
  theme(axis.title.y = element_text(face = "bold")) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#d94801") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
ot
#arrange in one figure
abun <- ggarrange(a,b,fi,fu,p,ot,
          ncol=3, nrow= 2)
abun
abun1 <- annotate_figure(abun, fig.lab = "C")

ggsave("output/plots/number_auxo_phylum_comparison.pdf", plot = abun1,
       width = 5, height = 4)
