##################   occurence of auxotrophies together ########################

Auxo <- Auxotrophy 
Auxo$Genomes <- NULL
Auxo$count <- NULL
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
View(occurence2)
occurence2 <- occurence[!(occurence$Var1 == 1 | occurence$Var2 == 1), ]
occurence2$Freq_all <-occurence2$Freq/nrow(Auxo)
occurence2$Freq_filt <- ifelse(occurence2$A1 == occurence2$A2, NA, occurence2$Freq/nrow(Auxo))

#occurence2$P <- occurence2$Freq/occurence2$perc

####################      visualization      ###################################
library(MetBrewer)

###Frequencies in relation to the number of genomes
 o <- ggplot(occurence2, aes(A1,A2, fill = Freq_filt)) +
  geom_tile(color ="white", lwd= 0.5, linetype = 1.5) +
  scale_fill_gradientn(colors = met.brewer("VanGogh3"), na.value = "gray94", limits = c(0, 0.4)) +
  guides(fill = guide_colourbar(barwidth = 15, barheight = 1, title ="Freq",
                                label = TRUE, ticks = FALSE)) +
  theme(legend.position = "bottom") +
  theme(axis.text.x = element_text(colour = "Black", size = 10, angle = 45, vjust = 0.5, hjust=0.5)) +
  theme(axis.text.y = element_text(colour = "black", size = 10, vjust = 0.5, hjust=0.5)) +
  theme(axis.title.x = element_text(colour = "Black", face = "bold", size = 10, margin = margin(5,0,0,0))) +
  theme(axis.title.y = element_text(colour = "Black", face = "bold", size = 10, margin = margin(0,5,0,0))) +
  theme(legend.title = element_text(colour = "black", size = 8, face = "bold", vjust = 1, margin = margin(0,10,1,0))) +
  theme(legend.text = element_text(size=8)) +
  xlab("Amino acid auxotrophy 1") +
  ylab("Amino acid auxotrophy 2") + 
  theme(panel.background = element_blank()) +
   scale_x_discrete(position = "top") +
  theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))
o
o1 <- annotate_figure(o, fig.lab = "B")

ggsave("output/plots/Heatmap_occurence_auxos.pdf", plot = o1,
       width = 5, height = 5)


# i <- ggplot(occurence2, aes(A1,A2, fill = P)) +
#   geom_tile(color ="white", lwd= 0.5, linetype = 1.5) +
#   scale_fill_gradientn(colors = met.brewer("OKeeffe2"), na.value = "gray94") +
#   guides(fill = guide_colourbar(barwidth = 15, barheight = 1, title ="Freq",
#                                 label = TRUE, ticks = FALSE)) +
#   theme(legend.position = "bottom") +
#   theme(axis.text.x = element_text(colour = "Black", size = 10, angle = 45, vjust = 0.5, hjust=0.5)) +
#   theme(axis.text.y = element_text(colour = "black", size = 10, vjust = 0.5, hjust=0.5)) +
#   theme(axis.title.x = element_text(colour = "Black", face = "bold", size = 12, margin = margin(5,0,0,0))) +
#   theme(axis.title.y = element_text(colour = "Black", face = "bold", size = 12, margin = margin(0,5,0,0))) +
#   theme(legend.title = element_text(colour = "black", size = 12, vjust = 1, margin = margin(0,10,1,0))) +
#   xlab("AA1") +
#   ylab("AA2") + 
#   theme(panel.background = element_blank()) +
#   theme(plot.margin= margin(0.5,0.5,0.5,0.5, "cm"))
# i





