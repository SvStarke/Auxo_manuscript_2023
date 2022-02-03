############# proportion of auxotrophic/prototrophic MAGs per phylum ###########

Auxotrophy_12 <- Auxotrophy_12[complete.cases(Auxotrophy_12), ]


a1 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 0 & Auxotrophy_12$phylum == "Actinobacteriota"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Actinobacteriota"])
a2 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 1 & Auxotrophy_12$phylum == "Actinobacteriota"])/ nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Actinobacteriota"])
a3 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 0 & Auxotrophy_12$phylum == "Bacteroidota"])/ nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Bacteroidota"])
a4 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 1 & Auxotrophy_12$phylum == "Bacteroidota"])/ nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Bacteroidota"])
a5 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 0 & Auxotrophy_12$phylum == "Other"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Other"]) 
a6 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 1 & Auxotrophy_12$phylum == "Other"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Other"])
a7 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 0 & Auxotrophy_12$phylum == "Firmicutes"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Firmicutes"])
a8 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 1 & Auxotrophy_12$phylum == "Firmicutes"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Firmicutes"])
a9 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 0 & Auxotrophy_12$phylum == "Fusobacteriota"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Fusobacteriota"])
a10 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 1 & Auxotrophy_12$phylum == "Fusobacteriota"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Fusobacteriota"])
a11 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 0 & Auxotrophy_12$phylum == "Proteobacteria"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Proteobacteria"])
a12 <- nrow(Auxotrophy_12[Auxotrophy_12$Status == 1 & Auxotrophy_12$phylum == "Proteobacteria"]) / nrow(Auxotrophy_12[Auxotrophy_12$phylum == "Proteobacteria"])

a <- data.table(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)
v <- t(a)
v1 <- data.table(v)
v1$Phylum <- c("Actinobacteriota", "Actinobacteriota", "Bacteroidota", "Bacteroidota",
               "Other",  "Other", "Firmicutes", "Firmicutes", "Fusobacteriota" , "Fusobacteriota", "Proteobacteria", "Proteobacteria")
v1$Status <- c("Auxotroph","Prototroph","Auxotroph","Prototroph","Auxotroph","Prototroph","Auxotroph","Prototroph","Auxotroph","Prototroph","Auxotroph","Prototroph")

####################         visualization         #############################

p <- ggplot(v1, aes(Phylum, V1)) +
  geom_bar(stat = "identity", aes(fill = Status)) +
  ylab("Proportion of MAGs") +
  scale_fill_manual(values = c('#fdd0a2','#fdae6b')) +
  theme(legend.position = "right") +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  theme(axis.text.x = element_text(size = 12, angle = 50, colour = "black", hjust = 1, margin = margin(0,0,0,0))) +
  theme(axis.text.y = element_text(size = 12, colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.ticks = element_blank())+
  theme(axis.title.y = element_text(colour = "Black", face = "bold", size = 14, margin = margin(0,10,0,0))) +
  theme(legend.title = element_text(colour = "black", size = 12, face = "bold",  vjust = 1, margin = margin(10,5,5,0))) +
  theme(legend.text = element_text(size = 12, colour = "black")) 
p

ggsave("/Users/Svenja/Desktop/Auxos_Proto_prop.pdf", plot = p,
       width = 7, height = 5)
  




