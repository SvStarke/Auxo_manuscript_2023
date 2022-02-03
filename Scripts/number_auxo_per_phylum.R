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
remove(numb)

#Firmicutes
fi <- ggplot(numb[phylum == "Firmicutes"], aes (count,abun, fill = "#fdae6b")) +
  geom_bar(stat = "identity") +
  xlab("Number of auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#fdae6b") +
  ggtitle("Firmicutes") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 8)) +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))
fi

#Actinobacteria
a <- ggplot(numb[phylum == "Actinobacteriota"], aes (count, abun, fill = "#8csd04")) +
  geom_bar(stat = "identity") +
  xlab("Number of auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") + 
  scale_fill_manual(values = "#8c2d04") +
  ggtitle("Actinobacteriota") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 8)) +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))
a

#Bacteroidetes
b <- ggplot(numb[phylum == "Bacteroidota"], aes (count,abun, fill = "#f16913")) +
  geom_bar(stat = "identity") +
  xlab("Number of auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#f16913") +
  ggtitle("Bacteroidota") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 8)) +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))
b

#Proteobacteria
p <- ggplot(numb[phylum == "Proteobacteria"], aes (count, abun, fill = "#fdd0a2")) +
  geom_bar(stat = "identity") +
  xlab("Number of auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = "#fdd0a2") +
  ggtitle("Proteobacteria") +
  coord_cartesian(ylim=c(0,60)) +
  theme(title = element_text(size = 8)) +
  theme(axis.text.x = element_text(colour = "black")) +
  theme(axis.text.y = element_text(colour = "black"))
p
#arrange in one figure
abun <- ggarrange(fi,a,b,p,
          labels = c("A","B","C", "D"),
          ncol=2, nrow= 2)

ggsave("/Users/Svenja/Desktop/number_auxo_phylum_comparison.pdf", plot = abun,
       width = 6, height = 5)
