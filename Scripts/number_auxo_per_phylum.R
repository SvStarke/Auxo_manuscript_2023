######################      numbers of auxotrophies per phylum #################

count <- unique(Auxotrophy_13$count)
phylum <- unique(Auxotrophy_13$phylum)
tables <- list()
k <- 1
for (i in count) {
  print(i) 
    for (pi in phylum) {
      x <- Auxotrophy_13[count == i & phylum == pi]
      table <- nrow(x)
      tablea <- data.frame(table)
      tablea$phylum <- pi
      tablea$count <- i
      tablea$countall <- nrow(Auxotrophy_13[Auxotrophy_13$phylum == pi])
      tables[[k]] <- tablea
      k <- k + 1
  }
}
numb <- rbindlist(tables)
numb$abun <- numb$table / numb$countall * 100

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

###number auxo per phylum and order
#Firmicutes
Auxotrophy_13_Firm <- Auxotrophy_13[Auxotrophy_13$phylum == "Firmicutes"]
relorder_Firm <- unique(names(sort(summary(as.factor(Auxotrophy_13_Firm$order)), decreasing = T)[1:10]))
rel_order_Firm <- unique(Auxotrophy_13$order)

#Bacteroidota
Auxotrophy_13_Bac <- Auxotrophy_13[Auxotrophy_13$phylum == "Bacteroidota"]
relorder_Bac <- names(sort(summary(as.factor(Auxotrophy_13_Bac$order)), decreasing = T)[1:10])
#remove NA values
relorder_Bac <- relorder_Bac[!is.na(relorder_Bac)]

#Actinobacteriota
Auxotrophy_13_Act <- Auxotrophy_13[Auxotrophy_13$phylum == "Actinobacteriota"]
relorder_Act <- names(sort(summary(as.factor(Auxotrophy_13_Act$order)), decreasing = T)[1:10])
#remove NA values
relorder_Act <- relorder_Act[!is.na(relorder_Act)]

#Fusobacteriota
Auxotrophy_13_Fus <- Auxotrophy_13[Auxotrophy_13$phylum == "Fusobacteriota"]
relorder_Fus <- names(sort(summary(as.factor(Auxotrophy_13_Fus$order)), decreasing = T)[1:10])
#remove NA values
relorder_Fus <- relorder_Fus[!is.na(relorder_Fus)]

#Proteobacteriota
Auxotrophy_13_Pro <- Auxotrophy_13[Auxotrophy_13$phylum == "Proteobacteria"]
relorder_Pro <- names(sort(summary(as.factor(Auxotrophy_13_Pro$order)), decreasing = T)[1:10])
#remove NA values
relorder_Pro <- relorder_Pro[!is.na(relorder_Pro)]


###
count <- unique(Auxotrophy_13_Firm$count)
Auxotrophy_13_Firm_filt <- subset(Auxotrophy_13_Firm, order %in% relorder_Firm)
Firm_list <- list()
k <- 1
for (i in count) {
  print(i) 
  for (oi in relorder_Firm) {
    x <- Auxotrophy_13_Firm[count == i & order == oi]
    nrows <- nrow(x)
    tax <- data.frame(nrows)
    tax$order <- oi
    tax$count <- i
    tax$countall <- nrow(Auxotrophy_13_Firm_filt[Auxotrophy_13_Firm_filt$count == i])
    Firm_list[[k]] <- tax
    k <- k + 1
  }
}

numb_Firm <- rbindlist(Firm_list)

rm(numb_Firm)
rm(Firm_list)
numb_Firm$abun <- numb_Firm$nrows / numb_Firm$countall * 100
View(numb_Firm)

#control calculated number of rows per order
sum(numb_Firm[which(numb_Firm$order == "Lachnospirales"), 1])

#visualization
fi_or <- ggplot(numb_Firm, aes (count,abun, fill = order)) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size=8)) +
  theme(legend.title = element_text(size=8)) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  coord_cartesian() +
  ggtitle("Firmicutes") +
  theme(title = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(axis.line.x = element_line(colour= "black"))+
  theme(axis.line.y = element_line(colour = "black"))+
  theme(title = element_text(size = 10)) +
  scale_fill_manual(values = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))

  
fi_or

####Bacteroidota
count <- unique(Auxotrophy_13_Bac$count)
Bac_list <- list()
k <- 1
for (i in count) {
  print(i) 
  for (oi in relorder_Bac) {
    x <- Auxotrophy_13_Bac[count == i & order == oi]
    nrows <- nrow(x)
    tax <- data.frame(nrows)
    tax$order <- oi
    tax$count <- i
    tax$countall <- nrow(Auxotrophy_13_Bac[Auxotrophy_13_Bac$count == i])
    Bac_list[[k]] <- tax
    k <- k + 1
  }
}

numb_Bac <- rbindlist(Bac_list)
rm(Bac_list)
numb_Bac$abun <- numb_Bac$nrows / numb_Bac$countall * 100

#visualization
ba_or <- ggplot(numb_Bac, aes (count,abun, fill = order)) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size=8)) +
  theme(legend.title = element_text(size=8)) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  coord_cartesian() +
  ggtitle("Bacteroidota") +
  theme(axis.line.x = element_line(colour= "black"))+
  theme(axis.line.y = element_line(colour = "black"))+
  theme(title = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(title = element_text(size = 10)) +
  scale_fill_manual(values = c("lightgrey","darkgrey")) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
ba_or



####Actinobacteriota
count <- unique(Auxotrophy_13_Act$count)
Act_list <- list()
k <- 1
for (i in count) {
  print(i) 
  for (oi in relorder_Act) {
    x <- Auxotrophy_13_Act[count == i & order == oi]
    nrows <- nrow(x)
    tax <- data.frame(nrows)
    tax$order <- oi
    tax$count <- i
    tax$countall <- nrow(Auxotrophy_13_Act[Auxotrophy_13_Act$count == i])
    Act_list[[k]] <- tax
    k <- k + 1
  }
}

numb_Act <- rbindlist(Act_list)
rm(Act_list)
numb_Act$abun <- numb_Act$nrows / numb_Act$countall * 100

#visualization
ac_or <- ggplot(numb_Act, aes (count,abun, fill = order)) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size=8)) +
  theme(legend.title = element_text(size=8)) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  coord_cartesian() +
  ggtitle("Actinobacteriota") +
  theme(axis.line.x = element_line(colour= "black"))+
  theme(axis.line.y = element_line(colour = "black"))+
  theme(title = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(title = element_text(size = 10)) +
  scale_fill_manual(values = c("#fee0b6","#d8daeb","#b2abd2","#8073ac","#542788")) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
ac_or



####Fusobacteriota
count <- unique(Auxotrophy_13_Fus$count)
Fus_list <- list()
k <- 1
for (i in count) {
  print(i) 
  for (oi in relorder_Fus) {
    x <- Auxotrophy_13_Fus[count == i & order == oi]
    nrows <- nrow(x)
    tax <- data.frame(nrows)
    tax$order <- oi
    tax$count <- i
    tax$countall <- nrow(Auxotrophy_13_Fus[Auxotrophy_13_Fus$count == i])
    Fus_list[[k]] <- tax
    k <- k + 1
  }
}

numb_Fus <- rbindlist(Fus_list)
rm(Fus_list)
numb_Fus$abun <- numb_Fus$nrows / numb_Fus$countall * 100

#visualization
fu_or <- ggplot(numb_Fus, aes (count,abun, fill = order)) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size=8)) +
  theme(legend.title = element_text(size=8)) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  coord_cartesian() +
  ggtitle("Fusobacteriota") +
  theme(axis.line.x = element_line(colour= "black"))+
  theme(axis.line.y = element_line(colour = "black"))+
  theme(title = element_text(size = 10)) +
  scale_fill_manual(values = c("black")) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
fu_or
#7f3b08
#b35806
#e08214
#fdb863

####Proteobacteriota
count <- unique(Auxotrophy_13_Pro$count)
Pro_list <- list()
k <- 1
for (i in count) {
  print(i) 
  for (oi in relorder_Pro) {
    x <- Auxotrophy_13_Pro[count == i & order == oi]
    nrows <- nrow(x)
    tax <- data.frame(nrows)
    tax$order <- oi
    tax$count <- i
    tax$countall <- nrow(Auxotrophy_13_Pro[Auxotrophy_13_Pro$count == i])
    Pro_list[[k]] <- tax
    k <- k + 1
  }
}

numb_Pro <- rbindlist(Pro_list)
rm(Pro_list)
numb_Pro$abun <- numb_Pro$nrows / numb_Pro$countall * 100

#visualization
pr_or <- ggplot(numb_Pro, aes (count,abun, fill = order)) +
  geom_bar(stat = "identity") +
  xlab("Auxotrophies per genome") +
  ylab("Abundance [%]") +
  theme(legend.position = "right") +
  theme(legend.text = element_text(size=8)) +
  theme(legend.title = element_text(size=8)) +
  theme(panel.background = element_rect(fill="white", colour= "white")) +
  coord_cartesian() +
  ggtitle("Proteobacteriota") +
  theme(axis.line.x = element_line(colour= "black"))+
  theme(title = element_text(face="bold")) +
  theme(axis.title.x = element_text(face="bold")) +
  theme(axis.title.y = element_text(face="bold")) +
  theme(axis.line.y = element_line(colour = "black"))+
  theme(title = element_text(size = 10)) +
  scale_fill_manual(values = c("#543005",
    "#8c510a",
    "#bf812d",
    "#dfc27d",
    "#f6e8c3",
    "#c7eae5",
    "#80cdc1",
    "#35978f",
    "#01665e",
    "#003c30")) +
  theme(axis.text.x = element_text(colour = "black", size = 10)) +
  theme(axis.text.y = element_text(colour = "black", size = 10))
pr_or

