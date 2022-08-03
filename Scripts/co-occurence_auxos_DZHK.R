########## co-occurence of auxos per person

Auxo <- Auxotrophy 


u <- merge(dzhk_relabun, Auxo, by.x="model", by.y="Genomes")
u$prop <- NULL


Auxo$Genomes <- NULL

samp <- unique(u$sample)
AA <- unique(colnames(Auxo))
AA2 <- unique(colnames(Auxo))
l <- list()
k <- 1
for(si in samp) {
  print(si)
    for(AAi1 in AA) {
      for(AAi2 in AA2) {
        sampi <- u[sample == si]
        nauxo <- table(sampi[[AAi1]],sampi[[AAi2]])
        occu <- data.frame(nauxo)
        occu$A1 <- AAi1
        occu$A2 <- AAi2
        occu$perc <- nrow(sampi[sampi[[AAi1]]== 0])
        occu <- data.table(occu)
        occu <- occu[Var1 == 0 & Var2 == 0]
        occu$numbrows <- nrow(sampi)
        occu$perc <- occu$Freq/occu$numbrows
        l[[k]] <- occu
        k <- k + 1
    }
     
  }
}

op <- rbindlist(l)
z <- aggregate(op$perc, by=list(AA1 = op$A1, AA2 =op$A2), FUN = sum)
z$perc <- z$x/190
z <- data.table(z)
z <- z[AA1 < AA2]
z <- z[AA2 != "Gly"]
z<- z[AA1 != "Gly"]
View(z)

oc <- ggplot(z, aes(AA1,AA2, fill = perc)) +
  geom_tile(color ="white", lwd= 0.5, linetype = 1.5) +
  scale_fill_gradientn(colors = met.brewer("VanGogh3"), limits = c(0, 0.6)) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1, title ="Freq",
                                label = TRUE, ticks = FALSE)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(colour = "Black", size = 10, angle = 45, vjust = 0.5, hjust=0.5)) +
  theme(axis.text.y = element_text(colour = "black", size = 10, vjust = 0.5, hjust=0.5)) +
  theme(axis.title.x = element_text(colour = "Black", size = 10, margin = margin(5,0,0,0))) +
  theme(axis.title.y = element_text(colour = "Black", size = 10, margin = margin(0,5,0,0))) +
  theme(legend.title = element_text(colour = "black", size = 8, vjust = 1, margin = margin(0,10,1,0))) +
  theme(legend.text = element_text(size=8)) +
  xlab("Amino acid auxotrophy 1") +
  ylab("Amino acid auxotrophy 2") + 
  theme(panel.background = element_blank()) +
  scale_x_discrete(position = "top") +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))
oc

ggsave("output/plots/Heatmap_occurence_auxos_persubject.pdf", plot = oc,
       width = 5, height = 5)

